# ============================================================================
# STANDALONE DE TEST FOR GSE296027 (run SEPARATELY from the pooled 5-dataset test)
# ============================================================================
# WHY THIS EXISTS:
# GSE296027 is FACS-sorted, plate-pooled microglia (~190 cells/sample,
# 63K-1.1M total UMI) vs. the other 5 datasets' droplet snRNA-seq of whole
# tissue (~3M-92M total UMI, 2-20K cells/sample). Pooling it into one TMM
# normalization with the other datasets makes calcNormFactors fail for all
# 7 of its samples (you'll see this as identical fallback norm factors and
# "no non-missing arguments to max" warnings in the pooled run). Including
# a ~GSE + group term does NOT fix this, because normalization factors are
# computed globally before the model is fit.
#
# This script re-normalizes GSE296027's 7 samples against ONLY each other,
# giving it a small but honest standalone AD-vs-Control test (n=4 vs n=3).
# Report this as a separate result, not merged into the pooled 5-dataset
# number from 05b_REAL_edgeR_analysis_FIXED.R.
# ============================================================================

if (!require("edgeR", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.rstudio.com")
  BiocManager::install("edgeR")
}
library(edgeR)

AD_PROJECT_DIR <- "/content/drive/MyDrive/AD_project"

counts_df <- read.csv(file.path(AD_PROJECT_DIR, "pseudobulk_counts_matrix.csv"),
                       row.names = 1, check.names = FALSE)
metadata  <- read.csv(file.path(AD_PROJECT_DIR, "pseudobulk_sample_metadata.csv"))

meta_296027 <- metadata[metadata$GSE == "GSE296027", ]
cat("GSE296027 samples:\n")
print(meta_296027[, c("sample", "condition", "n_cells", "total_UMI")])

counts_296027 <- as.matrix(counts_df[, meta_296027$sample])
# drop genes that are all-zero within this subset before filtering/normalizing
counts_296027 <- counts_296027[rowSums(counts_296027) > 0, ]
cat(paste("\nGenes with any signal in GSE296027 subset:", nrow(counts_296027), "\n"))

group <- factor(meta_296027$condition, levels = c("Control", "AD"))
cat("Groups:\n"); print(table(group))

y <- DGEList(counts = counts_296027, group = group)
cpm_mat <- cpm(y, log = FALSE)
keep <- rowSums(cpm_mat > 1) >= 2
y <- y[keep, , keep.lib.sizes = FALSE]
cat(paste("After filtering low-count genes:", nrow(y), "genes retained\n"))

y <- normLibSizes(y)
cat("Normalization factors (within-GSE296027 only):\n")
print(y$samples$norm.factors)

design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = ncol(design))

deg_296027 <- topTags(qlf, n = Inf)$table
names(deg_296027)[names(deg_296027) == "logFC"]  <- "log2FC"
names(deg_296027)[names(deg_296027) == "PValue"] <- "pvalue"
names(deg_296027)[names(deg_296027) == "FDR"]    <- "padj"

cat(paste("\nSignificant (padj<0.05):", sum(deg_296027$padj < 0.05, na.rm = TRUE),
          "of", nrow(deg_296027), "tested\n"))
cat("NOTE: with only 7 samples (4 AD, 3 Control), this test is underpowered.\n")
cat("Treat it as a standalone, exploratory result for this assay, not as\n")
cat("confirmatory evidence, and say so explicitly in the manuscript.\n")

out_file <- file.path(AD_PROJECT_DIR, "deg_results_GSE296027_standalone.csv")
write.csv(deg_296027, out_file, row.names = TRUE)
cat(paste("\n✓ Saved:", out_file, "\n"))
