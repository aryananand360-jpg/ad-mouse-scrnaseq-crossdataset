# ================================================================
# BATCH scRNA-seq CLUSTERING ANALYSIS
# 20 Mouse Alzheimer's Disease GEO Datasets
#
# Adapted from the single-GSE (GSE311935 HCL) reference pipeline.
# Runs the full workflow independently for each GSE ID below:
#
# GEO accession
#   -> Download supplementary files
#   -> Detect scRNA-seq expression matrix
#   -> QC (mouse mitochondrial genes: "mt-" prefix)
#   -> Normalization
#   -> Highly Variable Genes
#   -> PCA
#   -> Nearest-neighbor graph
#   -> Leiden clustering
#   -> UMAP
#   -> Marker genes (incl. mouse AD / microglia / DAM markers)
#   -> Save results
#
# Each GSE is processed inside a try/except block so that one
# failing dataset (e.g. no readable matrix, weird file layout)
# does not stop the rest of the batch.
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
import traceback
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
# 3. PROJECT SETTINGS
# ================================================================
GSE_IDS = [
    "GSE98969",
    "GSE140399",
    "GSE140510",
    "GSE143758",
    "GSE167600",
    "GSE178304",
    "GSE212606",
    "GSE218352",
    "GSE224398",
    "GSE224647",
    "GSE214244",
    "GSE233601",
    "GSE255662",
    "GSE285694",
    "GSE296027",
    "GSE164507",
    "GSE190012",
    "GSE264553",
    "GSE140511",
    "GSE98971",
]

ROOT_DIR = Path("/content/mouse_AD_scRNAseq")
ROOT_DIR.mkdir(parents=True, exist_ok=True)

# QC thresholds (adjust per-dataset if needed after inspecting violin plots)
MIN_GENES = 200
MAX_GENES = 6000
MAX_PCT_MT = 20

# Mouse Alzheimer's-relevant marker genes (mouse gene symbols use
# Title Case, e.g. "Cx3cr1", not "CX3CR1")
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

print("Datasets queued:", len(GSE_IDS))
print("Root directory:", ROOT_DIR)


# ================================================================
# HELPER FUNCTIONS
# ================================================================

def download_geo_supplementary_files(gse_id, suppl_dir):
    """Fetch the GEO suppl/ directory listing and download every file."""
    ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/suppl/"
    print("GEO supplementary directory:", ftp_base)

    response = requests.get(ftp_base, timeout=60)
    if response.status_code != 200:
        raise RuntimeError(
            f"Could not access GEO supplementary directory: {response.status_code}"
        )

    html = response.text
    files = re.findall(r'href="([^"]+)"', html)
    supp_files = [f for f in files if not f.endswith("/")]

    print("Files detected from GEO:")
    for f in supp_files:
        print(" ", f)

    downloaded_files = []
    for filename in supp_files:
        if filename in ["../", "?C=N;O=D"]:
            continue
        url = ftp_base + filename
        output_file = suppl_dir / filename
        print("Downloading:", filename)
        try:
            r = requests.get(url, stream=True, timeout=180)
            if r.status_code == 200:
                with open(output_file, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)
                downloaded_files.append(output_file)
                print(
                    "  Saved:", output_file,
                    f"({output_file.stat().st_size / 1024 / 1024:.2f} MB)"
                )
            else:
                print("  Download failed:", r.status_code)
        except Exception as e:
            print("  ERROR:", e)

    return downloaded_files


def extract_archives(suppl_dir):
    """Extract any tar/zip archives found in the supplementary directory."""
    print("\n---- Extracting archives ----")
    for f in list(suppl_dir.iterdir()):
        try:
            if f.is_file() and tarfile.is_tarfile(f):
                print("Extracting:", f.name)
                with tarfile.open(f) as tar:
                    tar.extractall(suppl_dir)
            elif f.is_file() and zipfile.is_zipfile(f):
                print("Extracting:", f.name)
                with zipfile.ZipFile(f) as z:
                    z.extractall(suppl_dir)
        except Exception as e:
            print("Could not extract:", f.name, e)


def find_expression_matrix(suppl_dir, base_dir):
    """Search for a 10x/H5/H5AD matrix and load it into an AnnData object."""
    matrix_files = list(suppl_dir.rglob("*matrix.mtx*"))
    h5_files = list(suppl_dir.rglob("*.h5"))
    h5ad_files = list(suppl_dir.rglob("*.h5ad"))
    rds_files = list(suppl_dir.rglob("*.rds"))

    print("MTX files:", [str(f) for f in matrix_files])
    print("H5 files:", [str(f) for f in h5_files])
    print("H5AD files:", [str(f) for f in h5ad_files])
    print("RDS files (not auto-loadable in Python):", [str(f) for f in rds_files])

    adata = None

    # OPTION A: H5AD
    if len(h5ad_files) > 0:
        print("Loading H5AD file...")
        adata = sc.read_h5ad(h5ad_files[0])

    # OPTION B: 10x H5
    elif len(h5_files) > 0:
        print("Loading 10x H5 file...")
        try:
            adata = sc.read_10x_h5(h5_files[0])
        except Exception as e:
            print("H5 loading failed:", e)

    # OPTION C: 10x MTX
    elif len(matrix_files) > 0:
        matrix_file = matrix_files[0]
        matrix_dir = matrix_file.parent
        print("10x MTX matrix detected in:", matrix_dir)

        barcode_candidates = (
            list(matrix_dir.glob("*barcodes.tsv*"))
            + list(matrix_dir.glob("*barcode*.tsv*"))
        )
        feature_candidates = (
            list(matrix_dir.glob("*features.tsv*"))
            + list(matrix_dir.glob("*genes.tsv*"))
            + list(matrix_dir.glob("*feature*.tsv*"))
        )

        print("Barcode files:", [str(f) for f in barcode_candidates])
        print("Feature files:", [str(f) for f in feature_candidates])

        if len(barcode_candidates) > 0 and len(feature_candidates) > 0:
            standard_dir = base_dir / "10x_matrix"
            standard_dir.mkdir(exist_ok=True)
            shutil.copy(matrix_file, standard_dir / "matrix.mtx.gz")
            shutil.copy(barcode_candidates[0], standard_dir / "barcodes.tsv.gz")
            shutil.copy(feature_candidates[0], standard_dir / "features.tsv.gz")

            print("Loading 10x matrix...")
            adata = sc.read_10x_mtx(
                standard_dir, var_names="gene_symbols", make_unique=True
            )

    if adata is None:
        print("\nNO DIRECT 10x/H5/H5AD EXPRESSION MATRIX WAS AUTOMATICALLY FOUND")
        print("Files present:")
        for f in suppl_dir.rglob("*"):
            if f.is_file():
                print(" ", f)
        raise RuntimeError(
            "No directly readable 10x/H5/H5AD matrix was found automatically. "
            "Inspect the printed files and add a custom loader for this GSE "
            "(e.g. a plain CSV/TSV counts table or an .rds file needing R)."
        )

    return adata


def process_gse(gse_id, root_dir):
    """Run the full clustering pipeline for a single mouse AD GSE dataset."""
    print("\n" + "=" * 70)
    print("PROCESSING:", gse_id)
    print("=" * 70)

    base_dir = root_dir / gse_id
    suppl_dir = base_dir / "supplementary"
    result_dir = base_dir / "results"
    suppl_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    # ---- GEO metadata ----
    print("\nDownloading GEO metadata...")
    gse = GEOparse.get_GEO(geo=gse_id, destdir=str(base_dir), silent=True)
    print("Number of samples:", len(gse.gsms))

    sample_table = []
    for gsm_id, gsm in gse.gsms.items():
        sample_table.append({
            "GSM": gsm_id,
            "Title": gsm.metadata.get("title", [""])[0],
            "Source": gsm.metadata.get("source_name_ch1", [""])[0],
        })
    sample_df = pd.DataFrame(sample_table)
    sample_df.to_csv(result_dir / f"{gse_id}_sample_info.csv", index=False)
    print(sample_df)

    # ---- Download + extract supplementary files ----
    print("\nDownloading supplementary files...")
    download_geo_supplementary_files(gse_id, suppl_dir)
    extract_archives(suppl_dir)

    # ---- Load expression matrix ----
    adata = find_expression_matrix(suppl_dir, base_dir)
    adata.var_names_make_unique()

    print("\nDataset info:")
    print(adata)
    print("Number of cells:", adata.n_obs)
    print("Number of genes:", adata.n_vars)

    # ---- QC (mouse mitochondrial genes are lowercase "mt-") ----
    print("\n---- Quality control ----")
    adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    qc_columns = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    print(adata.obs[qc_columns].describe())

    sc.pl.violin(
        adata, qc_columns, jitter=0.4, multi_panel=True,
        save=f"_{gse_id}_qc_violin.png", show=False
    )
    plt.close("all")

    sc.pl.scatter(
        adata, x="total_counts", y="n_genes_by_counts",
        save=f"_{gse_id}_counts_vs_genes.png", show=False
    )
    plt.close("all")

    sc.pl.scatter(
        adata, x="total_counts", y="pct_counts_mt",
        save=f"_{gse_id}_counts_vs_mt.png", show=False
    )
    plt.close("all")

    # Move scanpy's default figures/ output into this GSE's result dir
    for f in Path("figures").glob(f"*{gse_id}*"):
        shutil.move(str(f), result_dir / f.name)

    # ---- Filter low-quality cells ----
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
    print("Cells removed         :", before_cells - after_cells)

    if adata.n_obs < 50:
        raise RuntimeError(
            f"Only {adata.n_obs} cells survived QC for {gse_id} — "
            "too few for reliable clustering. Check matrix / QC thresholds."
        )

    # ---- Save raw counts, normalize ----
    adata.layers["counts"] = adata.X.copy()

    print("\nNormalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ---- Highly variable genes ----
    print("\nIdentifying highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    print("Highly variable genes:", adata.var["highly_variable"].sum())

    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()

    # ---- Scale + PCA ----
    print("\nScaling data and running PCA...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")

    # ---- Neighbors, Leiden, UMAP ----
    print("\nConstructing nearest-neighbor graph...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

    print("Running Leiden clustering...")
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden", random_state=42)
    print("Cluster sizes:")
    print(adata.obs["leiden"].value_counts().sort_index())

    print("Calculating UMAP...")
    sc.tl.umap(adata, random_state=42)

    sc.pl.umap(
        adata, color="leiden", legend_loc="on data", frameon=False,
        title=f"{gse_id} Mouse AD - Leiden Clusters",
        save=f"_{gse_id}_umap_leiden.png", show=False
    )
    plt.close("all")
    for f in Path("figures").glob(f"*{gse_id}*"):
        shutil.move(str(f), result_dir / f.name)

    # ---- Cluster composition ----
    cluster_counts = adata.obs["leiden"].value_counts().sort_index()
    cluster_percent = (
        adata.obs["leiden"].value_counts(normalize=True).sort_index() * 100
    )
    cluster_summary = pd.DataFrame({
        "Cell_Count": cluster_counts,
        "Percentage": cluster_percent.round(2),
    })
    cluster_summary.to_csv(result_dir / f"{gse_id}_cluster_composition.csv")

    # ---- Marker genes ----
    print("\nIdentifying marker genes...")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", pts=True)
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(result_dir / f"{gse_id}_all_marker_genes.csv", index=False)

    significant_markers = markers[
        (markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0.5)
    ].copy()
    top_markers = (
        significant_markers
        .sort_values(["group", "logfoldchanges"], ascending=[True, False])
        .groupby("group")
        .head(10)
    )
    top_markers.to_csv(
        result_dir / f"{gse_id}_top10_marker_genes_per_cluster.csv", index=False
    )
    print("Significant markers:", len(significant_markers))

    # ---- Mouse AD / microglia marker check + plots ----
    available_markers = [g for g in CANDIDATE_MARKERS if g in adata.raw.var_names]
    print("AD/microglia markers available in this dataset:", available_markers)

    if len(available_markers) > 0:
        sc.pl.umap(
            adata, color=available_markers, ncols=3, frameon=False,
            save=f"_{gse_id}_marker_umaps.png", show=False
        )
        plt.close("all")

        sc.pl.dotplot(
            adata, var_names=available_markers, groupby="leiden",
            standard_scale="var",
            save=f"_{gse_id}_marker_dotplot.png", show=False
        )
        plt.close("all")

        for f in Path("figures").glob(f"*{gse_id}*"):
            shutil.move(str(f), result_dir / f.name)

    # ---- Save AnnData + cell metadata + UMAP coords ----
    adata.write(base_dir / f"{gse_id}_clustered.h5ad")

    adata.obs.to_csv(result_dir / f"{gse_id}_cell_cluster_assignments.csv")

    umap_df = pd.DataFrame(
        adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP1", "UMAP2"]
    )
    umap_df["Cluster"] = adata.obs["leiden"].values
    umap_df.to_csv(result_dir / f"{gse_id}_UMAP_coordinates.csv")

    summary = {
        "gse_id": gse_id,
        "n_cells_final": adata.n_obs,
        "n_genes_used": adata.n_vars,
        "n_clusters": adata.obs["leiden"].nunique(),
        "status": "success",
    }
    print(f"\n{gse_id} complete —", summary)
    return summary


# ================================================================
# RUN THE BATCH
# ================================================================
batch_summary = []

for gse_id in GSE_IDS:
    try:
        result = process_gse(gse_id, ROOT_DIR)
        batch_summary.append(result)
    except Exception as e:
        print(f"\n*** {gse_id} FAILED: {e} ***")
        traceback.print_exc()
        batch_summary.append({
            "gse_id": gse_id,
            "n_cells_final": None,
            "n_genes_used": None,
            "n_clusters": None,
            "status": f"failed: {e}",
        })

# ================================================================
# FINAL BATCH SUMMARY
# ================================================================
summary_df = pd.DataFrame(batch_summary)
summary_path = ROOT_DIR / "batch_summary_mouse_AD_scRNAseq.csv"
summary_df.to_csv(summary_path, index=False)

print("\n" + "=" * 70)
print("BATCH RUN COMPLETE — 20 MOUSE ALZHEIMER'S DISEASE DATASETS")
print("=" * 70)
print(summary_df)
print("\nSummary saved to:", summary_path)
print("Per-dataset results are under:", ROOT_DIR, "/<GSE_ID>/results/")
