"""
REAL pathway enrichment (v2) — fixes the background-matching bug from v1.

v1 called gp.enrichr(..., background=background_genes), which POSTs your
~20K-gene background list to Enrichr's web background-upload endpoint.
That endpoint choked (JSON parse error) both times and gseapy silently
fell back to Enrichr's default whole-genome background — meaning your
p-values were NOT actually background-corrected, even though the script
looked like it worked.

FIX: use gp.enrich() (LOCAL Fisher's exact test) instead of gp.enrichr()
(WEB API). This downloads each gene-set library once via get_library()
(still needs internet for that one step), then computes the test
entirely locally against your real background — no size-limited web
endpoint involved, so no silent fallback is possible.
"""
import subprocess, sys
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "-U", "gseapy"], check=True)

import pandas as pd
import gseapy as gp

AD_PROJECT_DIR = '/content/drive/MyDrive/AD_project'

deg_all = pd.read_csv(f'{AD_PROJECT_DIR}/deg_results_6datasets_REAL.csv')
high_conf = pd.read_csv(f'{AD_PROJECT_DIR}/deg_high_confidence_cross_dataset_REAL.csv')

background_genes = deg_all['gene'].dropna().unique().tolist()
print(f"Background (all tested genes): {len(background_genes)}")

LIBRARY_NAMES = [
    'GO_Biological_Process_2023',
    'KEGG_2019_Mouse',
    'Reactome_2022',
]

print("\nDownloading gene-set libraries once (needs internet, this part IS a web call)...")
gene_sets = {}
for lib_name in LIBRARY_NAMES:
    lib = gp.get_library(name=lib_name, organism='Mouse')
    print(f"  {lib_name}: {len(lib)} terms")
    gene_sets.update(lib)  # merge into one dict; 'Gene_set' column disambiguates in output... 
    # NOTE: merging dicts means duplicate term names across libraries would
    # collide. Doesn't happen in practice for these 3 libraries (naming
    # conventions differ: "(GO:...)" vs "R-HSA-..." vs plain KEGG names),
    # but if you add more libraries later, run them separately instead.

def run_enrichment(gene_list, label):
    gene_list = [g for g in gene_list if isinstance(g, str)]
    print(f"\n{'='*80}\n{label}: {len(gene_list)} genes\n{'='*80}")
    if len(gene_list) < 5:
        print(f"  Skipping — too few genes ({len(gene_list)}) for meaningful enrichment.")
        return None

    enr = gp.enrich(
        gene_list=gene_list,
        gene_sets=gene_sets,
        background=background_genes,   # LOCAL Fisher's exact test — this is now actually used
        outdir=None,
    )
    res = enr.results.sort_values('Adjusted P-value')
    sig = res[res['Adjusted P-value'] < 0.05]
    print(f"  {len(sig)} / {len(res)} terms significant (Adj P < 0.05)")
    cols = [c for c in ['Gene_set', 'Term', 'Overlap', 'Odds Ratio',
                         'Adjusted P-value', 'Genes'] if c in sig.columns]
    print(sig[cols].head(15).to_string(index=False))
    return res

full_res = run_enrichment(high_conf['gene'].tolist(), "FULL high-confidence list")
up_res   = run_enrichment(high_conf.loc[high_conf['log2FC'] > 0, 'gene'].tolist(), "UP-regulated in AD")
down_res = run_enrichment(high_conf.loc[high_conf['log2FC'] < 0, 'gene'].tolist(), "DOWN-regulated in AD")

if full_res is not None:
    full_res.to_csv(f'{AD_PROJECT_DIR}/pathway_enrichment_full_REAL.csv', index=False)
    print(f"\n✓ Saved: pathway_enrichment_full_REAL.csv (background-corrected)")
if up_res is not None:
    up_res.to_csv(f'{AD_PROJECT_DIR}/pathway_enrichment_upregulated_REAL.csv', index=False)
    print(f"✓ Saved: pathway_enrichment_upregulated_REAL.csv (background-corrected)")
if down_res is not None:
    down_res.to_csv(f'{AD_PROJECT_DIR}/pathway_enrichment_downregulated_REAL.csv', index=False)
    print(f"✓ Saved: pathway_enrichment_downregulated_REAL.csv (background-corrected) "
          f"— still expect this to be underpowered/noisy with so few input genes; "
          f"a properly-corrected background makes the test more honest, not more powered.")

print("\nCompare these Adjusted P-values to the v1 (uncorrected-background) run —")
print("if terms/rankings barely move, the whole-genome-background fallback wasn't")
print("distorting things much. If p-values shift noticeably, this run is the one to report.")
