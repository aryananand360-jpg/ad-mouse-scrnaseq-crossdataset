"""
STEP 2 of 2 — combine per-dataset checkpoints into the final matrix.
======================================================================

Run this ONCE, after Step 1 has been run for every dataset you want
included (in however many separate sessions/runtimes that took).

This only reads the small per-dataset CSVs written by Step 1 — no
scanpy, no h5py, no .h5ad files — so it is cheap and safe to (re)run
any time, including in a completely fresh runtime, regardless of how
Step 1 was split across sessions. It will NOT crash the way the
combined single-script version did, because it never holds more than
one dataset's raw single-cell matrix in memory (it never touches raw
single-cell data at all).

If you re-run Step 1 for a dataset (e.g. to fix something), just re-run
this script afterward too — it always rebuilds from whatever is
currently in pseudobulk_parts/, so there's no stale-state risk.
"""

import glob
import os
import pandas as pd

AD_PROJECT_DIR = '/content/drive/MyDrive/AD_project'
PARTS_DIR = f'{AD_PROJECT_DIR}/pseudobulk_parts'

counts_files = sorted(glob.glob(f'{PARTS_DIR}/*_counts.csv'))
meta_files = sorted(glob.glob(f'{PARTS_DIR}/*_meta.csv'))

found_gses = sorted({os.path.basename(f).replace('_counts.csv', '') for f in counts_files})
print(f"Found checkpoints for {len(found_gses)} dataset(s): {found_gses}")

if len(found_gses) == 0:
    raise SystemExit(
        "No checkpoint files found in pseudobulk_parts/. "
        "Run 05a_step1_process_one_dataset.py for at least one dataset first."
    )

all_counts = {}   # gse -> DataFrame (genes x samples)
all_meta = []

for gse in found_gses:
    cdf = pd.read_csv(f'{PARTS_DIR}/{gse}_counts.csv', index_col=0)
    mdf = pd.read_csv(f'{PARTS_DIR}/{gse}_meta.csv')
    all_counts[gse] = cdf
    all_meta.append(mdf)
    print(f"  {gse}: {cdf.shape[1]} samples, {cdf.shape[0]} genes in its own panel")

# ---- align every dataset's genes onto the union of all gene symbols ----
union_genes = sorted(set().union(*[set(df.index) for df in all_counts.values()]))
print(f"\nUnion of gene symbols across {len(found_gses)} dataset(s): {len(union_genes)}")

aligned = []
for gse, df in all_counts.items():
    aligned.append(df.reindex(union_genes, fill_value=0))

counts_df = pd.concat(aligned, axis=1)
sample_df = pd.concat(all_meta, ignore_index=True)

# sanity: column order in counts_df must match sample_df['sample'] order
counts_df = counts_df[sample_df['sample'].tolist()]

print(f"\nFinal count matrix: {counts_df.shape} (genes x samples)")
print(f"\nCondition summary:\n{sample_df['condition'].value_counts()}")
print(f"\nCondition summary by dataset:")
print(pd.crosstab(sample_df['GSE'], sample_df['condition']))

# ---- flag genes not shared by every dataset (structural zeros from union-fill) ----
per_gse_gene_sets = {gse: set(df.index) for gse, df in all_counts.items()}
genes_in_all = set.intersection(*per_gse_gene_sets.values()) if len(per_gse_gene_sets) > 1 else set(union_genes)
print(f"\nGenes present in EVERY dataset's panel: {len(genes_in_all)} / {len(union_genes)}")
print("(Genes outside this set are 0 for any dataset that didn't measure them —")
print(" fine under the ~GSE + group design in 05b, but worth a second look if")
print(" a gene you care about isn't in this 'present in all' set.)")

partial_genes_file = f'{AD_PROJECT_DIR}/genes_not_in_all_datasets.csv'
pd.Series(sorted(set(union_genes) - genes_in_all), name='gene').to_csv(partial_genes_file, index=False)
print(f"  ✓ Saved list of partial-coverage genes: {partial_genes_file}")

# ---- save final combined files ----
counts_file = f'{AD_PROJECT_DIR}/pseudobulk_counts_matrix.csv'
meta_file = f'{AD_PROJECT_DIR}/pseudobulk_sample_metadata.csv'
counts_df.to_csv(counts_file)
sample_df.to_csv(meta_file, index=False)

print(f"\n✓ Saved: {counts_file}")
print(f"✓ Saved: {meta_file}")
print("\nReady for 05b_REAL_edgeR_analysis_FIXED.R")
