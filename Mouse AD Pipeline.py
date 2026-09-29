# ================================================================
# GSE_ID - Mouse Alzheimer's Disease
# Single-Cell RNA-seq Clustering Analysis
#
# Set GSE_ID below to any one of your 20 accessions and run the
# whole script. Re-run with a different GSE_ID for the next one.
#
# Workflow:
# GEO accession
# -> Download supplementary files
# -> Detect scRNA-seq expression matrix
# -> QC (mouse mitochondrial genes: "mt-" prefix)
# -> Normalization
# -> Highly Variable Genes
# -> PCA
# -> Nearest-neighbor graph
# -> Leiden clustering
# -> UMAP
# -> Marker genes (incl. mouse AD / microglia / DAM markers)
# -> Save results
# ================================================================

# ================================================================
# 1. INSTALL PACKAGES
# ================================================================
!pip -q install scanpy leidenalg igraph GEOparse requests beautifulsoup4

# ================================================================
# 2. IMPORT LIBRARIES
# ================================================================
import os
import re
import gzip
import shutil
import requests
import zipfile
import tarfile
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import GEOparse
from pathlib import Path

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
sc.settings.set_figure_params(figsize=(6, 5), dpi=100)
print("Scanpy version:", sc.__version__)

# ================================================================
# 3. PROJECT SETTINGS  <-- CHANGE GSE_ID HERE FOR EACH RUN
# ================================================================
GSE_ID = "GSE98969"

BASE_DIR = Path(f"/content/{GSE_ID}")
SUPPL_DIR = BASE_DIR / "supplementary"
RESULT_DIR = BASE_DIR / "results"
SUPPL_DIR.mkdir(parents=True, exist_ok=True)
RESULT_DIR.mkdir(parents=True, exist_ok=True)
print("Project directory:", BASE_DIR)

# QC thresholds (adjust after inspecting the violin plots if needed)
MIN_GENES = 200
MAX_GENES = 6000
MAX_PCT_MT = 20

# Mouse Alzheimer's-relevant marker genes (mouse symbols use Title Case,
# e.g. "Cx3cr1", not "CX3CR1")
CANDIDATE_MARKERS = [
    # Microglia / homeostatic
    "Cx3cr1", "P2ry12", "Tmem119", "Aif1", "Itgam", "Csf1r",
    # Disease-associated microglia (DAM) / AD-reactive microglia
    "Trem2", "Tyrobp", "Apoe", "Clu", "Cst7", "Lpl", "Itgax", "Cd9",
    "Axl", "Spp1",
    # Astrocytes
    "Gfap", "Aqp4", "Slc1a3", "Aldh1l1",
    # Neurons
    "Rbfox3", "Snap25", "Syt1", "Map2",
    # Oligodendrocytes / OPCs
    "Mbp", "Mog", "Plp1", "Pdgfra", "Olig1", "Olig2",
    # AD pathology genes
    "App", "Psen1", "Psen2", "Mapt", "Bin1",
]

# ================================================================
# 4. DOWNLOAD GEO SERIES INFORMATION
# ================================================================
print("\nDownloading GEO metadata...")
gse = GEOparse.get_GEO(geo=GSE_ID, destdir=str(BASE_DIR), silent=True)
print("GEO accession loaded:", GSE_ID)
print("\nNumber of samples:")
print(len(gse.gsms))

# ================================================================
# 5. DISPLAY SAMPLE INFORMATION
# ================================================================
print("\n================ SAMPLE INFORMATION ================\n")
sample_table = []
for gsm_id, gsm in gse.gsms.items():
    sample_table.append({
        "GSM": gsm_id,
        "Title": gsm.metadata.get("title", [""])[0],
        "Source": gsm.metadata.get("source_name_ch1", [""])[0]
    })
sample_df = pd.DataFrame(sample_table)
display(sample_df)
sample_df.to_csv(RESULT_DIR / f"{GSE_ID}_sample_info.csv", index=False)

# ================================================================
# 6. DOWNLOAD GEO SUPPLEMENTARY FILES
# ================================================================
print("\n================ DOWNLOADING SUPPLEMENTARY FILES ================\n")
ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GSE_ID[:-3]}nnn/{GSE_ID}/suppl/"
print("GEO supplementary directory:")
print(ftp_base)

response = requests.get(ftp_base)
if response.status_code != 200:
    raise RuntimeError(
        f"Could not access GEO supplementary directory: {response.status_code}"
    )

html = response.text
files = re.findall(r'href="([^"]+)"', html)
supp_files = [f for f in files if not f.endswith("/")]

print("\nFiles detected from GEO:")
for f in supp_files:
    print(" ", f)

# ================================================================
# 7. DOWNLOAD ALL SUPPLEMENTARY FILES
# ================================================================
print("\nDownloading files...\n")
downloaded_files = []
for filename in supp_files:
    if filename in ["../", "?C=N;O=D"]:
        continue
    url = ftp_base + filename
    output_file = SUPPL_DIR / filename
    print("Downloading:", filename)
    try:
        r = requests.get(url, stream=True, timeout=120)
        if r.status_code == 200:
            with open(output_file, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)
            downloaded_files.append(output_file)
            print(" Saved:", output_file,
                  f"({output_file.stat().st_size / 1024 / 1024:.2f} MB)")
        else:
            print(" Download failed:", r.status_code)
    except Exception as e:
        print(" ERROR:", e)

# ================================================================
# 8. LIST DOWNLOADED FILES
# ================================================================
print("\n================ DOWNLOADED FILES ================\n")
all_files = list(SUPPL_DIR.rglob("*"))
for f in all_files:
    if f.is_file():
        print(f.name, f"{f.stat().st_size / 1024 / 1024:.2f} MB")

# ================================================================
# 9. EXTRACT TAR / ZIP FILES IF PRESENT
# ================================================================
print("\n================ EXTRACTING ARCHIVES ================\n")
for f in list(SUPPL_DIR.iterdir()):
    try:
        if f.is_file() and tarfile.is_tarfile(f):
            print("Extracting:", f.name)
            with tarfile.open(f) as tar:
                tar.extractall(SUPPL_DIR)
        elif f.is_file() and zipfile.is_zipfile(f):
            print("Extracting:", f.name)
            with zipfile.ZipFile(f) as z:
                z.extractall(SUPPL_DIR)
    except Exception as e:
        print("Could not extract:", f.name, e)

# ================================================================
# 10. SEARCH FOR SINGLE-CELL EXPRESSION MATRICES
# ================================================================
print("\n================ SEARCHING FOR EXPRESSION MATRIX ================\n")
matrix_files = list(SUPPL_DIR.rglob("*matrix.mtx*"))
h5_files = list(SUPPL_DIR.rglob("*.h5"))
h5ad_files = list(SUPPL_DIR.rglob("*.h5ad"))
rds_files = list(SUPPL_DIR.rglob("*.rds"))

print("MTX files:")
for f in matrix_files:
    print(" ", f)
print("\nH5 files:")
for f in h5_files:
    print(" ", f)
print("\nH5AD files:")
for f in h5ad_files:
    print(" ", f)
print("\nRDS files (not auto-loadable in Python):")
for f in rds_files:
    print(" ", f)

# ================================================================
# 11. LOAD EXPRESSION DATA
# ================================================================
adata = None

# OPTION A: H5AD
if len(h5ad_files) > 0:
    print("\nLoading H5AD file...")
    adata = sc.read_h5ad(h5ad_files[0])

# OPTION B: 10x H5
elif len(h5_files) > 0:
    print("\nLoading 10x H5 file...")
    try:
        adata = sc.read_10x_h5(h5_files[0])
    except Exception as e:
        print("H5 loading failed:", e)

# OPTION C: 10x MTX
elif len(matrix_files) > 0:
    print("\n10x MTX matrix detected.")
    matrix_file = matrix_files[0]
    matrix_dir = matrix_file.parent
    print("Matrix directory:")
    print(matrix_dir)

    barcode_candidates = (
        list(matrix_dir.glob("*barcodes.tsv*"))
        + list(matrix_dir.glob("*barcode*.tsv*"))
    )
    feature_candidates = (
        list(matrix_dir.glob("*features.tsv*"))
        + list(matrix_dir.glob("*genes.tsv*"))
        + list(matrix_dir.glob("*feature*.tsv*"))
    )

    print("\nBarcode files:")
    for f in barcode_candidates:
        print(" ", f)
    print("\nFeature files:")
    for f in feature_candidates:
        print(" ", f)

    if len(barcode_candidates) > 0 and len(feature_candidates) > 0:
        standard_dir = BASE_DIR / "10x_matrix"
        standard_dir.mkdir(exist_ok=True)
        shutil.copy(matrix_file, standard_dir / "matrix.mtx.gz")
        shutil.copy(barcode_candidates[0], standard_dir / "barcodes.tsv.gz")
        shutil.copy(feature_candidates[0], standard_dir / "features.tsv.gz")

        print("\nLoading 10x matrix...")
        adata = sc.read_10x_mtx(
            standard_dir, var_names="gene_symbols", make_unique=True
        )

# IF NO MATRIX WAS FOUND
if adata is None:
    print("\n============================================================")
    print("NO DIRECT 10x/H5/H5AD EXPRESSION MATRIX WAS AUTOMATICALLY FOUND")
    print("============================================================")
    print("\nFiles downloaded from GEO:")
    for f in SUPPL_DIR.rglob("*"):
        if f.is_file():
            print(f)
    raise RuntimeError(
        "\nThe GEO supplementary files do not contain a directly "
        "readable 10x/H5/H5AD matrix using the automatic loader. "
        "Inspect the printed files and identify the expression matrix "
        "format before continuing (it may be a plain CSV/TSV counts "
        "table, or an .rds file that needs to be loaded in R)."
    )

# ================================================================
# 12. BASIC DATA INFORMATION
# ================================================================
print("\n================ DATASET INFORMATION ================\n")
print(adata)
print("\nNumber of cells:", adata.n_obs)
print("Number of genes:", adata.n_vars)

# ================================================================
# 13. MAKE GENE NAMES UNIQUE
# ================================================================
adata.var_names_make_unique()

# ================================================================
# 14. QUALITY CONTROL  (mouse mitochondrial genes: "mt-" prefix)
# ================================================================
print("\n================ QUALITY CONTROL ================\n")
adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

qc_columns = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
print(adata.obs[qc_columns].describe())

# ================================================================
# 15. QC VISUALIZATION
# ================================================================
sc.pl.violin(
    adata, qc_columns, jitter=0.4, multi_panel=True
)
plt.show()

# ================================================================
# 16. QC SCATTER PLOTS
# ================================================================
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
plt.show()

sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
plt.show()

# ================================================================
# 17. FILTER LOW QUALITY CELLS
# ================================================================
print("\nFiltering low-quality cells...")
before_cells = adata.n_obs
adata = adata[
    (adata.obs["n_genes_by_counts"] >= MIN_GENES)
    & (adata.obs["n_genes_by_counts"] <= MAX_GENES)
    & (adata.obs["pct_counts_mt"] < MAX_PCT_MT)
].copy()
after_cells = adata.n_obs
print("Cells before filtering:", before_cells)
print("Cells after filtering :", after_cells)
print("Cells removed :", before_cells - after_cells)

# ================================================================
# 18. SAVE RAW COUNTS
# ================================================================
adata.layers["counts"] = adata.X.copy()

# ================================================================
# 19. NORMALIZATION
# ================================================================
print("\nNormalizing data...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ================================================================
# 20. HIGHLY VARIABLE GENES
# ================================================================
print("\nIdentifying highly variable genes...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
print("Highly variable genes:", adata.var["highly_variable"].sum())

sc.pl.highly_variable_genes(adata)
plt.show()

# ================================================================
# 21. SAVE NORMALIZED FULL DATA
# ================================================================
adata.raw = adata

# ================================================================
# 22. KEEP HIGHLY VARIABLE GENES FOR PCA
# ================================================================
adata = adata[:, adata.var["highly_variable"]].copy()

# ================================================================
# 23. SCALE DATA
# ================================================================
print("\nScaling data...")
sc.pp.scale(adata, max_value=10)

# ================================================================
# 24. PCA
# ================================================================
print("\nRunning PCA...")
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")

# ================================================================
# 25. PCA VARIANCE PLOT
# ================================================================
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
plt.show()

# ================================================================
# 26. NEAREST-NEIGHBOR GRAPH
# ================================================================
print("\nConstructing nearest-neighbor graph...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# ================================================================
# 27. LEIDEN CLUSTERING
# ================================================================
print("\nRunning Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5, key_added="leiden", random_state=42)
print("\nCluster sizes:")
print(adata.obs["leiden"].value_counts().sort_index())

# ================================================================
# 28. UMAP
# ================================================================
print("\nCalculating UMAP...")
sc.tl.umap(adata, random_state=42)

# ================================================================
# 29. UMAP CLUSTER PLOT
# ================================================================
print("\n================ UMAP CLUSTERS ================\n")
sc.pl.umap(
    adata, color="leiden", legend_loc="on data", frameon=False,
    title=f"{GSE_ID} Mouse AD - Leiden Clusters"
)
plt.show()

# ================================================================
# 30. CLUSTER COMPOSITION
# ================================================================
cluster_counts = adata.obs["leiden"].value_counts().sort_index()
cluster_percent = adata.obs["leiden"].value_counts(normalize=True).sort_index() * 100

cluster_summary = pd.DataFrame({
    "Cell_Count": cluster_counts,
    "Percentage": cluster_percent.round(2)
})
print("\nCluster composition:")
display(cluster_summary)

# ================================================================
# 31. CLUSTER COMPOSITION PLOT
# ================================================================
cluster_percent.plot(kind="bar", figsize=(8, 5))
plt.xlabel("Leiden Cluster")
plt.ylabel("Percentage of Cells")
plt.title(f"{GSE_ID} Cluster Composition")
plt.tight_layout()
plt.show()

# ================================================================
# 32. FIND MARKER GENES
# ================================================================
print("\n================ IDENTIFYING MARKER GENES ================\n")
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", pts=True)

# ================================================================
# 33. TOP MARKER GENES
# ================================================================
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)
plt.show()

# ================================================================
# 34. EXTRACT MARKER GENES
# ================================================================
markers = sc.get.rank_genes_groups_df(adata, group=None)
print("\nTop marker genes:")
display(markers.head(30))

# ================================================================
# 35. FILTER SIGNIFICANT MARKERS
# ================================================================
significant_markers = markers[
    (markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0.5)
].copy()
print("\nNumber of significant markers:", len(significant_markers))

# ================================================================
# 36. TOP 10 MARKERS PER CLUSTER
# ================================================================
top_markers = (
    significant_markers
    .sort_values(["group", "logfoldchanges"], ascending=[True, False])
    .groupby("group")
    .head(10)
)
print("\n================ TOP MARKERS PER CLUSTER ================\n")
display(top_markers[["group", "names", "logfoldchanges", "pvals_adj"]])

# ================================================================
# 37. SAVE MARKER GENES
# ================================================================
markers.to_csv(RESULT_DIR / f"{GSE_ID}_all_marker_genes.csv", index=False)
top_markers.to_csv(
    RESULT_DIR / f"{GSE_ID}_top10_marker_genes_per_cluster.csv", index=False
)

# ================================================================
# 38. MOUSE ALZHEIMER'S / MICROGLIA MARKERS
# ================================================================
available_markers = [g for g in CANDIDATE_MARKERS if g in adata.raw.var_names]
print("\nMarkers available in dataset:")
print(available_markers)

# ================================================================
# 39. PLOT AVAILABLE BIOLOGICAL MARKERS
# ================================================================
if len(available_markers) > 0:
    sc.pl.umap(adata, color=available_markers, ncols=3, frameon=False)
    plt.show()

# ================================================================
# 40. DOTPLOT OF AD / MICROGLIA MARKERS
# ================================================================
if len(available_markers) > 0:
    sc.pl.dotplot(
        adata, var_names=available_markers, groupby="leiden",
        standard_scale="var"
    )
    plt.show()

# ================================================================
# 41. SAVE CLUSTERED ANNData OBJECT
# ================================================================
output_h5ad = BASE_DIR / f"{GSE_ID}_clustered.h5ad"
adata.write(output_h5ad)
print("\nAnnData saved:")
print(output_h5ad)

# ================================================================
# 42. SAVE CELL-LEVEL CLUSTER ASSIGNMENTS
# ================================================================
cell_metadata = adata.obs.copy()
cell_metadata.to_csv(RESULT_DIR / f"{GSE_ID}_cell_cluster_assignments.csv")

# ================================================================
# 43. SAVE UMAP COORDINATES
# ================================================================
umap_df = pd.DataFrame(
    adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP1", "UMAP2"]
)
umap_df["Cluster"] = adata.obs["leiden"].values
umap_df.to_csv(RESULT_DIR / f"{GSE_ID}_UMAP_coordinates.csv")

# ================================================================
# 44. FINAL SUMMARY
# ================================================================
print("\n")
print("=" * 70)
print(f"{GSE_ID} SINGLE-CELL RNA-seq ANALYSIS COMPLETED")
print("=" * 70)
print("\nDataset:")
print(" ", GSE_ID, "- Mouse Alzheimer's Disease")
print("\nCells:")
print(" ", adata.n_obs)
print("\nGenes used for clustering:")
print(" ", adata.n_vars)
print("\nLeiden clusters:")
print(" ", adata.obs["leiden"].nunique())
print("\nResults directory:")
print(RESULT_DIR)
print("\nGenerated files:")
for f in RESULT_DIR.iterdir():
    if f.is_file():
        print(" ", f.name)
print("\nMain output:")
print(f" {GSE_ID}_clustered.h5ad")
print("\nAnalysis complete.")
