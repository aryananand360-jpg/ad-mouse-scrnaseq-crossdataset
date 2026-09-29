#!/usr/bin/env Rscript

# ============================================================================
# REAL DIFFERENTIAL EXPRESSION TESTING WITH EDGER (FIXED)
# ============================================================================
#
# Fixes applied vs. the original 05b script:
#
# 1. BUG: `cat("="*80, "\n")` — R has no `*` operator for strings. This
#    line crashed the script before anything ran. Replaced with
#    strrep("=", 80).
# 2. BUG: `%+%` was used as a string-concatenation operator but is not a
#    base R operator and was never defined anywhere in the script — this
#    also crashed. Replaced with paste0().
# 3. DESIGN FIX: the original used exactTest(), which only supports a
#    single two-group factor. But your pseudobulk samples come from up
#    to 6 different datasets/technologies (see 05a output), so pooling
#    them into one AD-vs-Control test with no dataset term confounds
#    "dataset" with "disease status". This version builds
#    design = ~GSE + group  (glmQLFit/glmQLFTest) whenever there is more
#    than one dataset in the data, which controls for that. If you don't
#    have enough residual degrees of freedom for a full GSE term, the
#    script tells you and falls back to ~group with a printed warning
#    (do not silently trust that fallback in the manuscript text).
# 4. OPTION: set EXCLUDE_GSE below to drop datasets that use a
#    fundamentally different assay from the rest (e.g. GSE296027 is
#    plate-based Smart-seq of FACS-sorted microglia, vs. droplet
#    snRNA-seq of whole tissue in the others) as a sensitivity check.

# ---- CONFIG ----------------------------------------------------------------
EXCLUDE_GSE <- c()   # e.g. c("GSE296027") to run a tissue-only sensitivity model
# -----------------------------------------------------------------------------

sep80 <- function() cat(strrep("=", 80), "\n")

sep80()
cat("REAL DIFFERENTIAL EXPRESSION TESTING WITH EDGER (FIXED)\n")
sep80()
cat("\n")

if (!require("edgeR", quietly = TRUE)) {
  cat("Installing edgeR...\n")
  install.packages("BiocManager", repos = "http://cran.rstudio.com")
  BiocManager::install("edgeR")
}
library(edgeR)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("[1] Loading pseudobulk data...\n")

counts_file <- "/content/drive/MyDrive/AD_project/pseudobulk_counts_matrix.csv"
counts_df <- read.csv(counts_file, row.names = 1, check.names = FALSE)
cat(paste("  Count matrix:", nrow(counts_df), "genes x", ncol(counts_df), "samples\n"))

meta_file <- "/content/drive/MyDrive/AD_project/pseudobulk_sample_metadata.csv"
metadata <- read.csv(meta_file)
cat(paste("  Metadata:", nrow(metadata), "samples\n"))

if (length(EXCLUDE_GSE) > 0) {
  cat(paste("  Excluding datasets:", paste(EXCLUDE_GSE, collapse = ", "), "\n"))
  keep_samples <- !(metadata$GSE %in% EXCLUDE_GSE)
  metadata <- metadata[keep_samples, ]
  counts_df <- counts_df[, metadata$sample, drop = FALSE]
}

stopifnot(all(colnames(counts_df) == metadata$sample))

counts <- as.matrix(counts_df)

cat("\n  Samples per dataset x condition:\n")
print(table(metadata$GSE, metadata$condition))

# ============================================================================
# CREATE EDGER DGELIST
# ============================================================================

cat("\n[2] Creating DGEList object...\n")

group <- factor(metadata$condition, levels = c("Control", "AD"))
gse_factor <- factor(metadata$GSE)
cat("  Groups:\n")
print(table(group))

y <- DGEList(counts = counts, group = group)
cat("  DGEList created\n")

# ============================================================================
# QC & NORMALIZATION
# ============================================================================

cat("\n[3] QC and normalization...\n")

cpm_mat <- cpm(y, log = FALSE)
keep <- rowSums(cpm_mat > 1) >= 2
y <- y[keep, , keep.lib.sizes = FALSE]
cat(paste("  After filtering low-count genes:", nrow(y), "genes retained\n"))

y <- calcNormFactors(y)
cat("  TMM normalization complete\n")
cat("  Normalization factors:", round(y$samples$norm.factors, 3), "\n")

# ============================================================================
# DESIGN MATRIX — control for dataset-of-origin when >1 dataset present
# ============================================================================

cat("\n[4] Building design matrix...\n")

n_datasets <- nlevels(gse_factor)
use_gse_covariate <- FALSE

if (n_datasets > 1) {
  # Check every GSE has samples in more than one condition where possible,
  # and that we have enough residual df for a full ~GSE + group model.
  design_full <- tryCatch({
    model.matrix(~ gse_factor + group)
  }, error = function(e) NULL)

  df_resid <- ncol(y) - (if (!is.null(design_full)) ncol(design_full) else Inf)

  if (!is.null(design_full) && qr(design_full)$rank == ncol(design_full) && df_resid >= 1) {
    design <- design_full
    use_gse_covariate <- TRUE
    cat(paste("  Using design = ~GSE + group  (", n_datasets,
              "datasets,", df_resid, "residual df)\n"))
  } else {
    cat("  ⚠ WARNING: not enough residual degrees of freedom (or design is",
        "rank-deficient) to fit ~GSE + group with this sample size.\n")
    cat("  ⚠ Falling back to design = ~group. Dataset/batch effects are NOT",
        "controlled for in this run — treat results as preliminary and say",
        "so explicitly in the manuscript.\n")
    design <- model.matrix(~group)
  }
} else {
  cat("  Only one dataset present — using design = ~group\n")
  design <- model.matrix(~group)
}

print(design)

# ============================================================================
# DISPERSION ESTIMATION & TESTING (glmQLFit, not exactTest — supports covariates)
# ============================================================================

cat("\n[5] Estimating dispersion and fitting model...\n")

y <- estimateDisp(y, design)
cat(paste("  Common dispersion:", round(sqrt(y$common.dispersion), 4), "\n"))

fit <- glmQLFit(y, design)
coef_name <- tail(colnames(design), 1)  # the AD-vs-Control coefficient is always last
cat(paste("  Testing coefficient:", coef_name, "\n"))

qlf <- glmQLFTest(fit, coef = ncol(design))
cat("  ✓ Model fit and test complete\n")

deg_table <- topTags(qlf, n = Inf)$table
cat(paste("  Results for", nrow(deg_table), "genes\n"))

# ============================================================================
# POST-PROCESSING
# ============================================================================

cat("\n[6] Post-processing results...\n")

deg_table$gene <- rownames(deg_table)
names(deg_table)[names(deg_table) == "logFC"] <- "log2FC"
names(deg_table)[names(deg_table) == "PValue"] <- "pvalue"
names(deg_table)[names(deg_table) == "FDR"] <- "padj"

keep_cols <- intersect(c("gene", "log2FC", "logCPM", "pvalue", "padj"), names(deg_table))
deg_table <- deg_table[, keep_cols]
deg_table <- deg_table[order(deg_table$padj), ]

cat("✓ Results processed\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n")
sep80()
cat("RESULTS SUMMARY\n")
sep80()
cat("\n")

cat(paste("Model used:", if (use_gse_covariate) "~GSE + group (batch-controlled)" else "~group (NOT batch-controlled)", "\n"))
cat(paste("Total genes tested:", nrow(deg_table), "\n"))
cat(paste("Significant (padj < 0.05):", sum(deg_table$padj < 0.05, na.rm = TRUE), "\n"))
cat(paste("Significant (padj < 0.01):", sum(deg_table$padj < 0.01, na.rm = TRUE), "\n"))

sig <- deg_table$padj < 0.05
cat(paste("  Upregulated in AD:", sum(sig & deg_table$log2FC > 0, na.rm = TRUE), "\n"))
cat(paste("  Downregulated in AD:", sum(sig & deg_table$log2FC < 0, na.rm = TRUE), "\n"))

cat("\nTop 20 significant genes (ranked by padj):\n")
print(head(deg_table[which(deg_table$padj < 0.05), c("gene", "log2FC", "pvalue", "padj")], 20))

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n[7] Saving results...\n")

output_file <- "/content/drive/MyDrive/AD_project/deg_results_6datasets_REAL.csv"
write.csv(deg_table, output_file, row.names = FALSE)
cat(paste("  ✓ Full results saved to:", output_file, "\n"))

sig_file <- "/content/drive/MyDrive/AD_project/deg_significant_padj05_REAL.csv"
write.csv(deg_table[which(deg_table$padj < 0.05), ], sig_file, row.names = FALSE)
cat(paste("  ✓ Significant genes (padj<0.05) saved to:", sig_file, "\n"))

# ============================================================================
# QC PLOTS
# ============================================================================

cat("\n[8] Generating QC plots...\n")

pdf("/content/drive/MyDrive/AD_project/MA_plot_REAL.pdf", width = 8, height = 6)
plotMD(qlf, main = "MA Plot (AD vs Control, edgeR QL)")
abline(h = 0, col = "blue", lty = 2)
dev.off()
cat("  ✓ MA plot saved\n")

pdf("/content/drive/MyDrive/AD_project/volcano_plot_REAL.pdf", width = 8, height = 6)
with(deg_table, {
  plot(log2FC, -log10(pvalue),
       pch = 19, cex = 0.5, main = "Volcano Plot",
       xlab = "log2(FC) AD vs Control",
       ylab = "-log10(p-value)",
       col = ifelse(padj < 0.05, "red", "gray"))
  abline(v = c(-1, 1), h = -log10(0.05), lty = 2, col = "blue")
  legend("topright", legend = paste("Sig (padj<0.05):", sum(padj < 0.05, na.rm = TRUE)), bty = "n")
})
dev.off()
cat("  ✓ Volcano plot saved\n")

# ============================================================================
# VALIDATION & SANITY CHECKS
# ============================================================================

cat("\n[9] Validation checks...\n")

pval_valid <- all(deg_table$pvalue >= 0 & deg_table$pvalue <= 1, na.rm = TRUE)
padj_valid <- all(deg_table$padj >= 0 & deg_table$padj <= 1, na.rm = TRUE)
padj_geq_pval <- all(deg_table$padj >= deg_table$pvalue, na.rm = TRUE)
fc_numeric <- is.numeric(deg_table$log2FC)
not_all_sig <- mean(deg_table$padj < 0.05, na.rm = TRUE) < 0.5   # real data: not >50% "significant"
pvals_spread <- (max(deg_table$pvalue, na.rm=TRUE) - min(deg_table$pvalue, na.rm=TRUE)) > 0.5

cat(paste("  p-values valid (0-1 range):", pval_valid, "\n"))
cat(paste("  adjusted p-values valid (0-1 range):", padj_valid, "\n"))
cat(paste("  padj >= pvalue (as expected):", padj_geq_pval, "\n"))
cat(paste("  log2FC values are numeric:", fc_numeric, "\n"))
cat(paste("  fewer than 50% of genes 'significant':", not_all_sig, "\n"))
cat(paste("  p-values spread across a wide range:", pvals_spread, "\n"))

if (pval_valid & padj_valid & padj_geq_pval & fc_numeric & not_all_sig & pvals_spread) {
  cat("\n✓ ALL VALIDATION CHECKS PASSED — results look like real statistical output\n")
} else {
  cat("\n✗ ONE OR MORE VALIDATION CHECKS FAILED — investigate before trusting these numbers\n")
}

cat("\n")
sep80()
cat("END OF EDGER ANALYSIS\n")
sep80()
