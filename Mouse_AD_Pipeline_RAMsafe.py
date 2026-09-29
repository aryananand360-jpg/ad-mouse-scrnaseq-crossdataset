# ================================================================
# GSE_ID - Mouse Alzheimer's Disease
# MINIMAL, RAM-SAFE single-cell RNA-seq clustering pipeline
#
# What changed vs. the original script, and why:
#  - All plotting removed (plt.show() calls). Those figures already exist
#    in your original notebooks from a month ago; they cost RAM/time and
#    produce nothing new here.
#  - No permanent copy to Google Drive of downloads/extracted archives --
#    those stay on the fast local disk (/content/{GSE_ID}) exactly like
#    your original working run. Only the 3 files the integration steps
#    actually need get copied to Drive, at the very end.
#  - Clustering (scale + PCA + neighbors + leiden) runs on a TEMPORARY
#    copy restricted to the 2000 HVGs, which is deleted immediately after.
#    The full gene set + raw counts are kept in the main object the whole
#    time, because pseudobulk DE later needs ALL genes, not just the 2000
#    used for clustering.
#  - Prefers "filtered" (real-cell) matrices over "raw"/"unfiltered" ones
#    when GEO provides both -- loading the wrong one silently inflates
#    cell counts by 10-100x, which is the most likely reason QC/
#    normalization crashes before clustering even starts.
#  - gc.collect() + del after every stage that creates a large object.
#
# Set GSE_ID below and Runtime > Run all. Re-run with a different
# GSE_ID for each of your 9 datasets. Restart the runtime between
# datasets (Runtime > Restart runtime) so nothing lingers in memory.
# ================================================================

# ================================================================
# 1. INSTALL PACKAGES
# ================================================================
!pip -q install scanpy leidenalg igraph GEOparse requests

# ================================================================
# 2. IMPORTS
# ================================================================
import os
import re
import gc
import shutil
import requests
import zipfile
import tarfile
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import GEOparse
from pathlib import Path

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
print("Scanpy version:", sc.__version__)

# ================================================================
# 3. MOUNT DRIVE (only place Drive is touched -- final copy step)
# ================================================================
from google.colab import drive
drive.mount('/content/drive')

# ================================================================
# 4. PROJECT SETTINGS  <-- CHANGE GSE_ID HERE FOR EACH RUN
# ================================================================
GSE_ID = "GSE140510"

BASE_DIR = Path(f"/content/{GSE_ID}")          # fast local disk -- downloads, extraction, working files
SUPPL_DIR = BASE_DIR / "supplementary"
RESULT_DIR = BASE_DIR / "results"
SUPPL_DIR.mkdir(parents=True, exist_ok=True)
RESULT_DIR.mkdir(parents=True, exist_ok=True)

DRIVE_DIR = Path(f"/content/drive/MyDrive/AD_project/{GSE_ID}")  # only the 3 final files go here
DRIVE_DIR.mkdir(parents=True, exist_ok=True)

print("Local working directory:", BASE_DIR)
print("Drive output directory:", DRIVE_DIR)

MIN_GENES = 200
MAX_GENES = 6000
MAX_PCT_MT = 20
# Sanity-check threshold: if the loaded matrix has more "cells" than this,
# it's almost certainly an unfiltered/raw matrix (empty droplets included),
# not real cells. Stop before QC tries to process it and eats all your RAM.
MAX_PLAUSIBLE_CELLS = 300_000

# ================================================================
# 5. CHECK FOR AN EXISTING CHECKPOINT FIRST (resume, don't redo)
# ================================================================
qc_checkpoint = BASE_DIR / f"{GSE_ID}_qc_filtered.h5ad"
clustered_checkpoint = DRIVE_DIR / f"{GSE_ID}_clustered.h5ad"

if clustered_checkpoint.exists():
    print("This dataset is already fully done -- found:", clustered_checkpoint)
    print("Nothing to do. Skip to the next GSE_ID.")

# ================================================================
# 6. DOWNLOAD GEO SUPPLEMENTARY FILES (skip if already downloaded)
# ================================================================
if not qc_checkpoint.exists() and not clustered_checkpoint.exists():
    print("\nDownloading GEO metadata...")
    gse = GEOparse.get_GEO(geo=GSE_ID, destdir=str(BASE_DIR), silent=True)
    print("Number of samples:", len(gse.gsms))

    ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GSE_ID[:-3]}nnn/{GSE_ID}/suppl/"
    response = requests.get(ftp_base)
    if response.status_code != 200:
        raise RuntimeError(f"Could not access GEO supplementary directory: {response.status_code}")
    files = re.findall(r'href="([^"]+)"', response.text)
    supp_files = [f for f in files if not f.endswith("/") and f not in ("../", "?C=N;O=D")]

    for filename in supp_files:
        output_file = SUPPL_DIR / filename
        if output_file.exists():
            continue  # already downloaded in a previous attempt this session
        url = ftp_base + filename
        try:
            r = requests.get(url, stream=True, timeout=120)
            if r.status_code == 200:
                with open(output_file, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)
                print("Downloaded:", filename, f"({output_file.stat().st_size / 1024 / 1024:.2f} MB)")
            else:
                print("Download failed:", filename, r.status_code)
        except Exception as e:
            print("ERROR downloading", filename, ":", e)

    # ================================================================
    # 7. EXTRACT ARCHIVES
    # ================================================================
    for f in list(SUPPL_DIR.iterdir()):
        try:
            if f.is_file() and tarfile.is_tarfile(f):
                with tarfile.open(f) as tar:
                    tar.extractall(SUPPL_DIR)
            elif f.is_file() and zipfile.is_zipfile(f):
                with zipfile.ZipFile(f) as z:
                    z.extractall(SUPPL_DIR)
        except Exception as e:
            print("Could not extract:", f.name, e)

    # ================================================================
    # 8. FIND EXPRESSION MATRIX -- PREFER FILTERED OVER RAW/UNFILTERED
    # ================================================================
    all_matrix_files = list(SUPPL_DIR.rglob("*matrix.mtx*"))
    h5_files = list(SUPPL_DIR.rglob("*.h5"))
    h5ad_files = list(SUPPL_DIR.rglob("*.h5ad"))

    def is_raw_unfiltered(path):
        s = str(path).lower()
        return "raw" in s and "filtered" not in s

    filtered_only = [f for f in all_matrix_files if not is_raw_unfiltered(f)]
    matrix_files = filtered_only if filtered_only else all_matrix_files
    if all_matrix_files and not filtered_only:
        print("WARNING: only raw/unfiltered matrices found -- using them, but "
              "check cell count below carefully.")
    print(f"Using {len(matrix_files)} matrix file(s) "
          f"({'filtered' if filtered_only else 'raw'})")

    # ================================================================
    # 9. LOAD EXPRESSION DATA
    # ================================================================
    adata = None

    if len(h5ad_files) > 0:
        adata = sc.read_h5ad(h5ad_files[0])
    elif len(h5_files) > 0:
        adata = sc.read_10x_h5(h5_files[0])
    elif len(matrix_files) > 0:
        adatas = []
        for matrix_file in matrix_files:
            matrix_dir = matrix_file.parent
            sample_id = matrix_file.name.split('_matrix.mtx')[0]
            barcode_candidates = list(matrix_dir.glob(f"{sample_id}*barcodes.tsv*"))
            feature_candidates = (list(matrix_dir.glob(f"{sample_id}*features.tsv*")) +
                                  list(matrix_dir.glob(f"{sample_id}*genes.tsv*")))
            if barcode_candidates and feature_candidates:
                standard_sample_dir = BASE_DIR / f"10x_matrix/{sample_id}"
                standard_sample_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy(matrix_file, standard_sample_dir / "matrix.mtx.gz")
                shutil.copy(barcode_candidates[0], standard_sample_dir / "barcodes.tsv.gz")
                shutil.copy(feature_candidates[0], standard_sample_dir / "features.tsv.gz")
                sample_adata = sc.read_10x_mtx(standard_sample_dir, var_names="gene_symbols", make_unique=True)
                sample_adata.obs['sample_id'] = sample_id
                adatas.append(sample_adata)
        if adatas:
            adata = sc.concat(adatas, join="outer", uns_merge='first', label="sample_id",
                               keys=[a.obs['sample_id'].iloc[0] for a in adatas]) if len(adatas) > 1 else adatas[0]
            del adatas
            gc.collect()

    if adata is None:
        raise RuntimeError("No readable 10x/H5/H5AD matrix found. Inspect SUPPL_DIR manually.")

    # ================================================================
    # 10. SANITY CHECK BEFORE DOING ANYTHING EXPENSIVE
    # ================================================================
    print("Loaded shape (cells x genes):", adata.shape)
    if adata.n_obs > MAX_PLAUSIBLE_CELLS:
        raise RuntimeError(
            f"Loaded {adata.n_obs} 'cells' -- this is almost certainly a raw/"
            f"unfiltered matrix (empty droplets included), not real cells. "
            f"Stopping before QC to avoid a RAM crash. Inspect {SUPPL_DIR} and "
            f"identify the filtered/real-cell files manually."
        )

    adata.var_names_make_unique()

    # ================================================================
    # 11. QC AND FILTER
    # ================================================================
    adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    before_cells = adata.n_obs
    adata = adata[
        (adata.obs["n_genes_by_counts"] >= MIN_GENES)
        & (adata.obs["n_genes_by_counts"] <= MAX_GENES)
        & (adata.obs["pct_counts_mt"] < MAX_PCT_MT)
    ].copy()
    print("Cells before filtering:", before_cells, "| after:", adata.n_obs)

    adata.layers["counts"] = adata.X.copy()

    # Checkpoint right here -- if clustering crashes below, you resume from
    # this QC'd object instead of re-downloading and re-QC'ing from scratch.
    adata.write(qc_checkpoint)
    print("Checkpoint saved locally:", qc_checkpoint)

    gc.collect()

# ================================================================
# 12. RESUME POINT: load the QC checkpoint if we didn't just make it
# ================================================================
if 'adata' not in dir() or adata is None:
    adata = sc.read_h5ad(qc_checkpoint)
    print("Resumed from QC checkpoint:", qc_checkpoint, "| shape:", adata.shape)

# ================================================================
# 13. NORMALIZE + HVG (sparse, in place -- cheap, no densification yet)
# ================================================================
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
print("Highly variable genes:", adata.var["highly_variable"].sum())

# ================================================================
# 14. CLUSTERING ON A TEMPORARY HVG-ONLY COPY (deleted right after)
#     Full adata (all genes, raw counts layer) is untouched by this.
# ================================================================
tmp = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(tmp, max_value=10)
sc.tl.pca(tmp, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(tmp, n_neighbors=15, n_pcs=30)
sc.tl.leiden(tmp, resolution=0.5, key_added="leiden", random_state=42)
sc.tl.umap(tmp, random_state=42)

adata.obs["leiden"] = tmp.obs["leiden"].values
adata.obsm["X_pca"] = tmp.obsm["X_pca"]
adata.obsm["X_umap"] = tmp.obsm["X_umap"]

del tmp
gc.collect()
print("Cluster sizes:")
print(adata.obs["leiden"].value_counts().sort_index())

# ================================================================
# 15. MARKER GENES (on full gene set, log-normalized, not scaled)
# ================================================================
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=False, pts=True)
markers = sc.get.rank_genes_groups_df(adata, group=None)
significant_markers = markers[(markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0.5)].copy()
top_markers = (
    significant_markers.sort_values(["group", "logfoldchanges"], ascending=[True, False])
    .groupby("group").head(10)
)

# ================================================================
# 16. SAVE THE 3 FILES THE INTEGRATION STEPS ACTUALLY NEED
# ================================================================
output_h5ad = BASE_DIR / f"{GSE_ID}_clustered.h5ad"
adata.write(output_h5ad)

cell_metadata_path = RESULT_DIR / f"{GSE_ID}_cell_cluster_assignments.csv"
adata.obs.to_csv(cell_metadata_path)

markers_path = RESULT_DIR / f"{GSE_ID}_top10_marker_genes_per_cluster.csv"
top_markers.to_csv(markers_path, index=False)

# Copy just these 3 to Drive -- everything else stays local/disposable
shutil.copy(output_h5ad, DRIVE_DIR / output_h5ad.name)
shutil.copy(cell_metadata_path, DRIVE_DIR / cell_metadata_path.name)
shutil.copy(markers_path, DRIVE_DIR / markers_path.name)

print("\nDone. Saved to Drive:")
for f in DRIVE_DIR.iterdir():
    print(" ", f.name)
print(f"\nCells: {adata.n_obs} | Genes: {adata.n_vars} | Clusters: {adata.obs['leiden'].nunique()}")
