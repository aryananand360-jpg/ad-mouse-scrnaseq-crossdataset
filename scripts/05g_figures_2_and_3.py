"""
Figure 2 (single-cell landscape) and Figure 3 (cell-type composition, AD vs Control).

RUN ORDER: run fix_cell_annotation_confidence.py FIRST -- this script reads its
output (cell_annotation_master_6datasets_CONFIDENCE_TIERED.csv), not the raw
CORRECTED file, so zero-evidence labels are already demoted.

Low-RAM: reads only obs columns + the 2-D UMAP via h5py. Never loads counts or
builds an AnnData, so it is safe for GSE214244 (100k+ cells).

GSE296027 is left out of both figures: it is FACS-sorted microglia, so cell-type
proportions are meaningless there, and its AD/Control labels are plate-confounded.
"""
import os, re
import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AD_DIR = os.environ.get("AD_DIR", "/content/drive/MyDrive/AD_project")
OUT_DIR = f"{AD_DIR}/figures_REAL"
os.makedirs(OUT_DIR, exist_ok=True)
matplotlib.rcParams['pdf.fonttype'] = 42

DATASETS = ["GSE140399", "GSE140510", "GSE143758", "GSE208683", "GSE214244"]
NATIVE = {"GSE140399", "GSE140510", "GSE214244"}     # have obs['sample_id']
CROSSWALK = {                                        # barcode-prefix datasets
    "GSE143758": {
        'EZ.Batch1_HipR-AD-G1-4w': 'AD', 'EZ.Batch1_HipR-AD-G3-2w': 'AD',
        'EZ.Batch1_HipR-WT-G1-4w': 'Control', 'EZ.Batch1_HipR-WT-G3-2w': 'Control',
        'NP40.Batch1_Wt-Hip-S2-R': 'Control', 'NP40.Batch1_Untreated-Hip-S2-R': 'AD',
        'NP40.Batch3_Wt-Hip-S1-L': 'Control', 'NP40.Batch3_Wt-Hip-S2-L': 'Control',
        'NP40.Batch3_Untreated-Hip-S1-L': 'AD', 'NP40.Batch3_Untreated-Hip-S2-L': 'AD'},
    "GSE208683": {'WT_1': 'Control', 'WT_2': 'Control', 'AD_1': 'AD', 'AD_2': 'AD'},
}
MAX_POINTS = 40000   # UMAP downsample per dataset (plotting only)

# ---------------- h5py readers (same logic as 05a_step1) ----------------
def _dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])

def read_field(node):
    if isinstance(node, h5py.Group):
        k = set(node.keys())
        if {'categories', 'codes'} <= k:
            cats, codes = _dec(node['categories'][:]), node['codes'][:]
            out = np.empty(len(codes), dtype=object)
            ok = codes >= 0
            out[ok] = cats[codes[ok]]; out[~ok] = None
            return out
        if {'values', 'mask'} <= k:
            v = node['values'][:]
            if v.dtype.kind in ('S', 'O'): v = _dec(v)
            v = np.array(v, dtype=object)
            m = node['mask'][:].astype(bool)
            if m.mean() > 0.99:      # broken mask convention -> ignore it
                return v
            v[m] = None
            return v
        raise ValueError(f"unknown encoding {sorted(k)}")
    v = node[:]
    return _dec(v) if v.dtype.kind in ('S', 'O') else v

def read_obs(f, col):
    return read_field(f['obs'][col])

def read_index(f):
    name = f['obs'].attrs.get('_index', '_index')
    return read_field(f['obs'][name])

# ---------------- inputs ----------------
ann = pd.read_csv(f"{AD_DIR}/cell_annotation_master_6datasets_CONFIDENCE_TIERED.csv")
master = pd.read_csv(f"{AD_DIR}/master_metadata.csv")
lookup = {(str(r.GSE), str(r.GSM)): r.condition_FINAL for r in master.itertuples()}
GSM_RE = re.compile(r'(GSM\d+)')

PALETTE = {
    'Excitatory_neurons': '#d62728', 'Inhibitory_neurons': '#8c564b', 'Neurons': '#e377c2',
    'Neurons_Unclassified': '#f4a6a6', 'Astrocytes': '#1f77b4', 'Oligodendrocytes': '#bcbd22',
    'OPCs': '#7f7f7f', 'Endothelial': '#2ca02c', 'Homeostatic_microglia': '#9467bd',
    'Activated_DAM_microglia': '#17becf', 'Perivascular_macrophages': '#ff7f0e',
    'Choroid_plexus_ependymal': '#c49c94', 'Low_quality_high_mito': '#333333',
    'Ambiguous_needs_manual_review': '#cccccc',
}
def color_for(ct, _extra={}):
    if ct in PALETTE: return PALETTE[ct]
    if ct not in _extra:
        _extra[ct] = plt.cm.tab20(len(_extra) % 20)
    return _extra[ct]

def load_dataset(gse):
    """Returns DataFrame: cell_type, sample, condition (+ umap arrays)."""
    path = f"{AD_DIR}/{gse}/{gse}_clustered.h5ad"
    with h5py.File(path, 'r') as f:
        leiden = read_obs(f, 'leiden').astype(str)
        umap = f['obsm/X_umap'][:] if 'obsm' in f and 'X_umap' in f['obsm'] else None
        if gse in NATIVE:
            sid = read_obs(f, 'sample_id').astype(str)
            m = [GSM_RE.search(s) for s in sid]
            cond = np.array([lookup.get((gse, x.group(1)), 'Unknown') if x else 'Unknown' for x in m])
            sample = sid
        else:
            bcs = read_index(f)
            cw = CROSSWALK[gse]
            sample = np.array([next((p for p in cw if b is not None and b.startswith(p)), None)
                               for b in bcs], dtype=object)
            cond = np.array([cw.get(s, 'Unknown') if s else 'Unknown' for s in sample])
    a = ann[ann.GSE == gse]
    cmap = dict(zip(a.cluster.astype(str), a.final_cell_type))
    ct = np.array([cmap.get(l, 'Ambiguous_needs_manual_review') for l in leiden])
    df = pd.DataFrame({'cell_type': ct, 'sample': sample, 'condition': cond})
    return df, umap

# ---------------- run ----------------
data = {}
for g in DATASETS:
    if not os.path.exists(f"{AD_DIR}/{g}/{g}_clustered.h5ad"):
        print(f"skip {g}: file missing"); continue
    data[g] = load_dataset(g)
    print(f"{g}: {len(data[g][0])} cells")

# ===== FIGURE 2: UMAP landscape =====
fig, axes = plt.subplots(2, 3, figsize=(14, 8.5))
axes = axes.ravel()
seen = set()
rng = np.random.default_rng(0)
for ax, (g, (df, umap)) in zip(axes, data.items()):
    if umap is None:
        ax.text(0.5, 0.5, f"{g}\nno UMAP stored", ha='center'); ax.axis('off'); continue
    idx = np.arange(len(df))
    if len(idx) > MAX_POINTS: idx = rng.choice(idx, MAX_POINTS, replace=False)
    cols = [color_for(c) for c in df.cell_type.values[idx]]
    ax.scatter(umap[idx, 0], umap[idx, 1], s=0.6, c=cols, linewidths=0, rasterized=True)
    ax.set_title(f"{g} (n={len(df):,} cells)", fontsize=10)
    ax.axis('off')
    seen |= set(df.cell_type.unique())
legend_ax = axes[len(data)] if len(data) < 6 else None
if legend_ax is not None:
    legend_ax.axis('off')
    handles = [plt.Line2D([0], [0], marker='o', ls='', color=color_for(c), label=c) for c in sorted(seen)]
    legend_ax.legend(handles=handles, loc='center', fontsize=8, frameon=False, title="Cell type")
fig.suptitle("Figure 2. Single-cell landscape and cluster-level cell-type annotation", fontsize=13)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/Figure2_UMAP_landscape.pdf", dpi=300)
fig.savefig(f"{OUT_DIR}/Figure2_UMAP_landscape.png", dpi=300)
plt.close(fig)
print("saved Figure 2")

# ===== FIGURE 3: composition per sample =====
rows = []
for g, (df, _) in data.items():
    d = df[df.condition.isin(['AD', 'Control']) & df['sample'].notna()]
    ct = pd.crosstab(d['sample'], d['cell_type'])
    prop = ct.div(ct.sum(axis=1), axis=0)
    cond = d.groupby('sample')['condition'].first()
    for s in prop.index:
        for c in prop.columns:
            rows.append({'GSE': g, 'sample': s, 'condition': cond[s], 'cell_type': c,
                         'proportion': prop.loc[s, c], 'n_cells': int(ct.loc[s].sum())})
comp = pd.DataFrame(rows)
comp.to_csv(f"{AD_DIR}/cell_type_proportions_per_sample_REAL.csv", index=False)

# Show cell types making up >=1% of cells on average; the rest are in the CSV
mean_prop = comp.groupby('cell_type')['proportion'].mean().sort_values(ascending=False)
shown = mean_prop[mean_prop >= 0.01].index.tolist()
ncol = 4; nrow = int(np.ceil(len(shown) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.4 * nrow), squeeze=False)
gse_cols = dict(zip(DATASETS, plt.cm.Set1.colors))
for ax, ct in zip(axes.ravel(), shown):
    sub = comp[comp.cell_type == ct]
    for i, cond in enumerate(['Control', 'AD']):
        y = sub[sub.condition == cond]
        jitter = np.random.default_rng(1).uniform(-0.12, 0.12, len(y))
        ax.scatter(np.full(len(y), i) + jitter, y.proportion * 100,
                   c=[gse_cols[g] for g in y.GSE], s=28, edgecolor='k', linewidth=0.3)
        if len(y):
            ax.hlines(y.proportion.mean() * 100, i - 0.25, i + 0.25, color='k', lw=1.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['Control', 'AD'])
    ax.set_title(ct.replace('_', ' '), fontsize=9)
    ax.set_ylabel('% of cells in sample', fontsize=8)
for ax in axes.ravel()[len(shown):]:
    ax.axis('off')
handles = [plt.Line2D([0], [0], marker='o', ls='', color=gse_cols[g], label=g,
                      markeredgecolor='k') for g in data]
fig.legend(handles=handles, loc='lower right', fontsize=8, title="Dataset")
fig.suptitle("Figure 3. Cell-type composition by sample (each dot = one sample; bar = mean)\n"
             "Exploratory: n = 2-5 samples per group per dataset; no significance testing", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(f"{OUT_DIR}/Figure3_cell_type_composition.pdf", dpi=300)
fig.savefig(f"{OUT_DIR}/Figure3_cell_type_composition.png", dpi=300)
print("saved Figure 3; proportions CSV written")
