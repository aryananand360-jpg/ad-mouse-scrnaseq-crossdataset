"""
REAL cross-dataset consistency check — replaces the fabricated
cross_dataset_consistency_6datasets.csv (44 fake genes).

For each significant gene from the pooled model, checks whether the
AD-vs-Control direction is consistent when computed independently
within EACH dataset (using that dataset's own samples only, simple
mean log2(CPM+1) difference — not a formal test, just a directional
sanity check). This is what tells you whether "290 upregulated genes"
is a broad, real, cross-dataset signal, or an artifact of one dataset
with an outsized effect dominating the pooled ~GSE + group test.

GSE296027 is excluded here (see 05c_GSE296027_standalone.R for why) —
this checks consistency across the 5 comparable droplet-snRNA-seq
datasets only.
"""
import pandas as pd
import numpy as np

AD_PROJECT_DIR = '/content/drive/MyDrive/AD_project'

deg = pd.read_csv(f'{AD_PROJECT_DIR}/deg_results_6datasets_REAL.csv')
counts = pd.read_csv(f'{AD_PROJECT_DIR}/pseudobulk_counts_matrix.csv', index_col=0)
meta = pd.read_csv(f'{AD_PROJECT_DIR}/pseudobulk_sample_metadata.csv')
meta = meta[meta['GSE'] != 'GSE296027'].copy()  # excluded — see 05c script

sig_genes = deg.loc[deg['padj'] < 0.05, 'gene'].tolist()
print(f"Checking cross-dataset consistency for {len(sig_genes)} significant genes "
      f"across {meta['GSE'].nunique()} datasets ({meta['GSE'].unique().tolist()})...")

# library-size-normalize per sample (simple CPM) within the 5-dataset subset
lib_sizes = counts[meta['sample']].sum(axis=0)
cpm = counts[meta['sample']].div(lib_sizes, axis=1) * 1e6
log2cpm = np.log2(cpm + 1)

rows = []
for gene in sig_genes:
    if gene not in log2cpm.index:
        continue
    pooled_dir = np.sign(deg.loc[deg['gene'] == gene, 'log2FC'].iloc[0])

    per_gse_dir = []
    for gse, sub in meta.groupby('GSE'):
        ad_samples = sub.loc[sub['condition'] == 'AD', 'sample']
        ctrl_samples = sub.loc[sub['condition'] == 'Control', 'sample']
        if len(ad_samples) == 0 or len(ctrl_samples) == 0:
            continue
        diff = log2cpm.loc[gene, ad_samples].mean() - log2cpm.loc[gene, ctrl_samples].mean()
        per_gse_dir.append(np.sign(diff))

    n_tested = len(per_gse_dir)
    n_agree = sum(d == pooled_dir for d in per_gse_dir)
    rows.append({
        'gene': gene,
        'pooled_log2FC': deg.loc[deg['gene'] == gene, 'log2FC'].iloc[0],
        'pooled_padj': deg.loc[deg['gene'] == gene, 'padj'].iloc[0],
        'n_datasets_tested': n_tested,
        'n_datasets_same_direction': n_agree,
        'consensus_direction': 'AD_up' if pooled_dir > 0 else 'AD_down',
        'fraction_consistent': n_agree / n_tested if n_tested else np.nan,
    })

result = pd.DataFrame(rows).sort_values('fraction_consistent', ascending=False)
out_file = f'{AD_PROJECT_DIR}/cross_dataset_consistency_REAL.csv'
result.to_csv(out_file, index=False)

print(f"\nFraction of significant genes with 100% cross-dataset agreement: "
      f"{(result['fraction_consistent'] == 1.0).mean():.1%}")
print(f"Fraction with 5/5 datasets agreeing: "
      f"{((result['n_datasets_same_direction'] == 5) & (result['n_datasets_tested'] == 5)).mean():.1%}")
print(f"\nTop 15 by consistency:")
print(result.head(15).to_string(index=False))
print(f"\n✓ Saved: {out_file}")
print("\nIf most of your 'significant' genes only agree in 1-2 of 5 datasets,")
print("the pooled result is being driven by one dataset's effect size, not a")
print("broad real signal — investigate which dataset before reporting this.")
