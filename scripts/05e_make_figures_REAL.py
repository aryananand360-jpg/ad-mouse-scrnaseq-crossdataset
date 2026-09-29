"""
Figures 4, 5, 7 — built entirely from the REAL files produced so far.
(Figure 3, the volcano/MA plots, already exist from 05b_REAL_edgeR_analysis_FIXED.R)

Requires: pseudobulk_counts_matrix.csv, pseudobulk_sample_metadata.csv,
deg_results_6datasets_REAL.csv, deg_high_confidence_cross_dataset_REAL.csv,
pathway_enrichment_full_REAL.csv (all already on Drive from earlier steps).
"""
import subprocess, sys
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "matplotlib", "seaborn"], check=True)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # editable text in Illustrator/Inkscape
matplotlib.rcParams['font.size'] = 9

AD_PROJECT_DIR = '/content/drive/MyDrive/AD_project'
FIG_DIR = f'{AD_PROJECT_DIR}/figures_REAL'
import os
os.makedirs(FIG_DIR, exist_ok=True)

deg_all = pd.read_csv(f'{AD_PROJECT_DIR}/deg_results_6datasets_REAL.csv')
high_conf = pd.read_csv(f'{AD_PROJECT_DIR}/deg_high_confidence_cross_dataset_REAL.csv')
counts = pd.read_csv(f'{AD_PROJECT_DIR}/pseudobulk_counts_matrix.csv', index_col=0)
meta = pd.read_csv(f'{AD_PROJECT_DIR}/pseudobulk_sample_metadata.csv')
meta_5 = meta[meta['GSE'] != 'GSE296027'].copy()  # excluded — plate-confounded, see earlier notes

# ============================================================================
# FIGURE 4 — heatmap of high-confidence DEGs across the 5 comparable datasets
# ============================================================================
lib_sizes = counts[meta_5['sample']].sum(axis=0)
cpm = counts[meta_5['sample']].div(lib_sizes, axis=1) * 1e6
log2cpm = np.log2(cpm + 1)

top_genes = high_conf.sort_values('padj').head(50)['gene'].tolist()
top_genes = [g for g in top_genes if g in log2cpm.index]
mat = log2cpm.loc[top_genes]
z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0)

col_order = meta_5.sort_values(['condition', 'GSE'])['sample'].tolist()
z = z[col_order]
meta_ordered = meta_5.set_index('sample').loc[col_order]

fig, (ax_ann, ax_hm) = plt.subplots(
    2, 1, figsize=(12, 11), gridspec_kw={'height_ratios': [0.3, 10]}, sharex=True)

cond_colors = meta_ordered['condition'].map({'AD': '#c0392b', 'Control': '#2980b9'})
ax_ann.bar(range(len(col_order)), [1]*len(col_order), color=cond_colors, width=1.0)
ax_ann.set_xlim(-0.5, len(col_order)-0.5)
ax_ann.set_yticks([]); ax_ann.set_xticks([])
for spine in ax_ann.spines.values():
    spine.set_visible(False)
ax_ann.set_title('Figure 4. High-confidence DEGs (top 50 by padj), pseudobulk log2CPM, z-scored\n'
                  'Top bar: red = AD, blue = Control (5 datasets, GSE296027 excluded)', fontsize=10)

im = ax_hm.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
ax_hm.set_yticks(range(len(top_genes)))
ax_hm.set_yticklabels(top_genes, fontsize=6)
ax_hm.set_xticks(range(len(col_order)))
ax_hm.set_xticklabels(col_order, rotation=90, fontsize=5)
fig.colorbar(im, ax=ax_hm, label='z-score (log2CPM)', shrink=0.6)
plt.tight_layout()
fig.savefig(f'{FIG_DIR}/Figure4_DEG_heatmap_REAL.pdf', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"✓ Figure 4 saved: {FIG_DIR}/Figure4_DEG_heatmap_REAL.pdf")

# ============================================================================
# FIGURE 5 — DAM / homeostatic marker panel, per-dataset log2FC dot plot
# ============================================================================
dam_up = ['Trem2', 'Apoe', 'Lpl', 'Spp1', 'Cst7', 'Itgax', 'Axl', 'Ccl3', 'Ccl4', 'Lilrb4a']
homeostatic_down = ['P2ry12', 'Cx3cr1', 'Tmem119', 'Siglech', 'Selplg']
marker_genes = [g for g in dam_up + homeostatic_down if g in log2cpm.index]

rows = []
for gene in marker_genes:
    for gse, sub in meta_5.groupby('GSE'):
        ad_s = sub.loc[sub['condition'] == 'AD', 'sample']
        ct_s = sub.loc[sub['condition'] == 'Control', 'sample']
        if len(ad_s) == 0 or len(ct_s) == 0:
            continue
        l2fc = log2cpm.loc[gene, ad_s].mean() - log2cpm.loc[gene, ct_s].mean()
        rows.append({'gene': gene, 'GSE': gse, 'log2FC': l2fc,
                      'category': 'DAM (expect up)' if gene in dam_up else 'Homeostatic (expect down)'})
marker_df = pd.DataFrame(rows)
gene_order = dam_up + homeostatic_down
gene_order = [g for g in gene_order if g in marker_df['gene'].unique()]
gse_order = sorted(marker_df['GSE'].unique())

fig, ax = plt.subplots(figsize=(7, max(4, 0.4*len(gene_order))))
for i, gene in enumerate(gene_order):
    sub = marker_df[marker_df['gene'] == gene]
    for j, gse in enumerate(gse_order):
        row = sub[sub['GSE'] == gse]
        if len(row) == 0:
            continue
        val = row['log2FC'].iloc[0]
        color = '#c0392b' if val > 0 else '#2980b9'
        ax.scatter(j, i, s=abs(val)*80 + 20, color=color, alpha=0.75, edgecolor='k', linewidth=0.3)

ax.set_yticks(range(len(gene_order)))
ax.set_yticklabels(gene_order)
ax.set_xticks(range(len(gse_order)))
ax.set_xticklabels(gse_order, rotation=45, ha='right')
ax.axhline(len(dam_up) - 0.5, color='gray', linestyle='--', linewidth=0.8)
ax.set_title('Figure 5. DAM / homeostatic marker log2FC per dataset\n'
              'Red = up in AD, Blue = down in AD, size = |log2FC|', fontsize=10)
ax.set_xlim(-0.5, len(gse_order)-0.5)
plt.tight_layout()
fig.savefig(f'{FIG_DIR}/Figure5_marker_panel_REAL.pdf', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"✓ Figure 5 saved: {FIG_DIR}/Figure5_marker_panel_REAL.pdf")

# ============================================================================
# FIGURE 7 — pathway enrichment bar plot (top 15 terms, full high-confidence list)
# ============================================================================
enrich = pd.read_csv(f'{AD_PROJECT_DIR}/pathway_enrichment_full_REAL.csv')
top_terms = enrich.sort_values('Adjusted P-value').head(15).copy()
top_terms['neglog10padj'] = -np.log10(top_terms['Adjusted P-value'])
top_terms = top_terms.sort_values('neglog10padj')

palette = {'GO_Biological_Process_2023': '#27ae60', 'KEGG_2019_Mouse': '#e67e22',
           'Reactome_2022': '#8e44ad'}
colors = top_terms['Gene_set'].map(palette).fillna('#7f8c8d')

fig, ax = plt.subplots(figsize=(9, 6))
ax.barh(top_terms['Term'], top_terms['neglog10padj'], color=colors)
ax.set_xlabel('-log10(adjusted p-value)')
ax.set_title('Figure 7. Top enriched pathways — high-confidence DEGs (n=235)\n'
              'Background-corrected Fisher\'s exact test vs. all tested genes', fontsize=10)
handles = [plt.Rectangle((0,0),1,1, color=c) for c in palette.values()]
ax.legend(handles, palette.keys(), loc='lower right', fontsize=7)
plt.tight_layout()
fig.savefig(f'{FIG_DIR}/Figure7_pathway_enrichment_REAL.pdf', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"✓ Figure 7 saved: {FIG_DIR}/Figure7_pathway_enrichment_REAL.pdf")

print(f"\nAll figures saved to: {FIG_DIR}")
