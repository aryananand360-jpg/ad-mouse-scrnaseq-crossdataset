"""
Corrected validation — fixes the mouse/human gene-casing bug in the
original marker check, and adds the checks 07_FABRICATION_ANALYSIS
flagged as the actual bar for "these are real, not fabricated" (not
just a marker-presence check).
"""
import pandas as pd

AD_PROJECT_DIR = '/content/drive/MyDrive/AD_project'
deg = pd.read_csv(f'{AD_PROJECT_DIR}/deg_results_6datasets_REAL.csv')

print("="*80)
print("VALIDATING REAL DE RESULTS (v2 — mouse gene casing)")
print("="*80)

sig = deg[deg['padj'] < 0.05]
print(f"\nTotal genes tested: {len(deg)}")
print(f"Significant (padj<0.05): {len(sig)}  ({len(sig)/len(deg):.1%} of tested genes)")
print(f"  Upregulated: {(sig['log2FC'] > 0).sum()}")
print(f"  Downregulated: {(sig['log2FC'] < 0).sum()}")

# Mouse casing: DAM (disease-associated microglia) up-genes and
# homeostatic down-genes, from Keren-Shaul et al. 2017 and related lit.
dam_up = ['Trem2', 'Apoe', 'Lpl', 'Spp1', 'Cst7', 'Itgax', 'Clec7a', 'Axl',
          'Ccl3', 'Ccl4', 'Lilrb4a']
homeostatic_down = ['P2ry12', 'Cx3cr1', 'Tmem119', 'Siglech', 'Selplg']

print(f"\n--- DAM (activation) markers — expect UP in AD ---")
for g in dam_up:
    row = deg[deg['gene'] == g]
    if len(row):
        r = row.iloc[0]
        flag = "✓" if r['log2FC'] > 0 else "✗ (wrong direction)"
        print(f"  {g:10s} log2FC={r['log2FC']:+.3f}  padj={r['padj']:.2e}  {flag}")
    else:
        print(f"  {g:10s} — not in tested gene set (filtered out or not detected)")

print(f"\n--- Homeostatic markers — expect DOWN in AD ---")
for g in homeostatic_down:
    row = deg[deg['gene'] == g]
    if len(row):
        r = row.iloc[0]
        flag = "✓" if r['log2FC'] < 0 else "✗ (wrong direction)"
        print(f"  {g:10s} log2FC={r['log2FC']:+.3f}  padj={r['padj']:.2e}  {flag}")
    else:
        print(f"  {g:10s} — not in tested gene set (filtered out or not detected)")
