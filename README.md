# ad-mouse-scrnaseq-crossdataset
Cross-dataset reanalysis of mouse AD scRNA-seq data

# Reproducible Cross-Dataset Single-Cell Transcriptomic Signature of Alzheimer's Disease
A meta-analysis of six independent mouse brain scRNA-seq datasets (360,620 cells) 
using standardized preprocessing, pseudobulk differential expression, and cross-dataset 
consistency analysis.

## Datasets
- GSE140399, GSE140510, GSE143758, GSE208683, GSE214244, GSE296027

## Key Findings
- 301 significant differentially expressed genes (FDR < 0.05)
- 290 upregulated (96.3%), dominated by immune/interferon genes
- 235 genes (78%) show consistent direction across datasets
- Strong enrichment for Immune System (p = 1.9 × 10⁻⁴³) and Innate Immune pathways

## Pipeline
1. Quality control & clustering (Scanpy)
2. Confidence-tiered cell-type annotation (marker-based)
3. Sample-level pseudobulk aggregation
4. Differential expression (edgeR)
5. Cross-dataset consistency analysis
6. Pathway enrichment (Reactome, KEGG, GO)
7. Leave-one-dataset-out robustness validation

## Files
- `/scripts/` — Analysis pipelines and validation scripts
- `/results/` — Processed DEG tables, pseudobulk matrices, consistency/pathway results
- `/data/` — Cell-type annotations, supplementary tables

## Usage
```bash
# Example: run full pipeline for one dataset
python scripts/Mouse_AD_Pipeline_RAMsafe.py  # Set GSE_ID inside script

# Run validation checks
Rscript scripts/03_LEAVE_ONE_OUT_ROBUSTNESS.R
python scripts/04_LOOCV_VISUALIZATION.py
```

## Requirements
- Python 3.8+: scanpy, pandas, numpy, scipy
- R 4.0+: edgeR, tidyverse
- See scripts for full package list and versions

## Citation
[Authors, Year]. [Manuscript title]. *bioRxiv* [DOI/link].

## License
MIT
