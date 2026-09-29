"""
Builds the two gene lists that actually belong in the manuscript,
replacing the fabricated deg_results_6datasets.csv /
cross_dataset_consistency_6datasets.csv entirely.
"""
import pandas as pd

AD_PROJECT_DIR = '/content/drive/MyDrive/AD_project'

deg = pd.read_csv(f'{AD_PROJECT_DIR}/deg_results_6datasets_REAL.csv')
consistency = pd.read_csv(f'{AD_PROJECT_DIR}/cross_dataset_consistency_REAL.csv')

merged = deg.merge(consistency[['gene', 'n_datasets_tested', 'n_datasets_same_direction',
                                 'fraction_consistent']], on='gene', how='left')

sig = merged[merged['padj'] < 0.05].copy()
high_confidence = sig[sig['fraction_consistent'] == 1.0].copy()

print(f"Results 3.4 number — significant (padj<0.05), pooled 5-dataset model: {len(sig)}")
print(f"Results 3.5 number — of those, consistent direction in ALL 5 datasets: "
      f"{len(high_confidence)} ({len(high_confidence)/len(sig):.1%})")

sig = sig.sort_values('padj')
high_confidence = high_confidence.sort_values('padj')

sig.to_csv(f'{AD_PROJECT_DIR}/deg_significant_padj05_REAL_annotated.csv', index=False)
high_confidence.to_csv(f'{AD_PROJECT_DIR}/deg_high_confidence_cross_dataset_REAL.csv', index=False)

print(f"\n✓ Saved: deg_significant_padj05_REAL_annotated.csv  ({len(sig)} genes — use for Fig 3/4)")
print(f"✓ Saved: deg_high_confidence_cross_dataset_REAL.csv  ({len(high_confidence)} genes — use for Fig 5/7 and as your 'reproducible genes' number)")
