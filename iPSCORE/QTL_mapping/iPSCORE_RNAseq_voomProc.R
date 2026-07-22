library(edgeR) 
library(limma)
library(tidyverse)
library(data.table)

counts <- read.delim("~/ppc_gene_counts.txt", row.names = 1)

# Need to remove leading X from any column names
colnames(counts) <- gsub("X", "", colnames(counts))
# Need to replace "." with "-"
colnames(counts) <- gsub("\\.", "-", colnames(counts))

# Sample Metadata
count_meta <- read_csv("~/data/RNA_SampleMetaData.csv")

#colnames(counts) %in% count_meta$RNA_UUID # colnames of counts are RNA_UUID
# There are only 220 iPSC samples, the rest are either pancreatic progenitor cell samples
# or cardiovascular progenitor cell samples

# What are potential grouping variables?
# clone and passage, family ID

# Build a crosswalk BEFORE renaming columns
# Do this right after loading count_meta_iPSC

# Need to create a map from colnames to their conditions
count_meta_PPC <- count_meta %>% 
  filter(Tissue == "PPC")

id_crosswalk <- count_meta_PPC %>%
  filter(RNA_UUID %in% colnames(counts)) %>%   # only samples in your expression data
  select(RNA_UUID, WGS_UUID, Subject_UUID, iPSCORE_ID, Family_ID, Clone, Passage) %>%
  mutate(formatted_id = paste0("Fam", Family_ID, "_Clone", Clone, "_Pass", Passage))

# Save it for later
write_tsv(id_crosswalk, "~/id_crosswalk.tsv")

# See which Fam/Clone/Pass combos are duplicated
#dupes <- id_crosswalk$formatted_id[duplicated(id_crosswalk$formatted_id)]
#id_crosswalk %>% filter(formatted_id %in% dupes) %>% arrange(formatted_id)

# IT says that there are duplciates, but I'm not sure if that's actually the case
# THey do have identical family_id, clone, and passage, but not iPSCORE_ID

# head(count_meta_iPSC)
# 
# head(counts)

# Find the index of each colname in count_meta_iPSC$RNA_UUID
match_indices <- match(colnames(counts), count_meta_PPC$RNA_UUID)

# Create new column names using the matched indices
colnames(counts) <- paste0(
  "Fam",count_meta_PPC$Family_ID[match_indices],"_",
  "Clone", count_meta_PPC$Clone[match_indices],"_",
  "Pass",count_meta_PPC$Passage[match_indices]
)

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

# Split the column names to extract the three variables
split_names <- strsplit(colnames(counts), "_")

# Extract each component
# family <- sapply(split_names, function(x) x[1])
clone <- sapply(split_names, function(x) x[2])
passage <- sapply(split_names, function(x) x[3])

# Create the group vector using interaction
group <- interaction(clone, passage, drop = TRUE)

plotMDS(d, col = as.numeric(group))
# There's a clear line down the middle, but its not clear 

mm <- model.matrix(~0 + group)

y <- voom(d, mm, plot = T)

express <- y$E
exp_var <- 1/y$weights


## ============================================================
## GTEx-style Gene Expression Normalization & Filtering Pipeline
## ============================================================
## Steps:
##   1. TMM normalization (edgeR)
##   2. Autosomal gene filtering (≥0.1 TPM in ≥20% samples &
##      ≥6 raw reads in ≥20% samples)
##   3. Inverse normal transformation across samples per gene
## ============================================================

library(edgeR)

# -----------------------------------------------------------------
# 0. Load data
# -----------------------------------------------------------------
#counts <- read.delim("/work/jb621/GEMMAmod_iPSCORE/ppc_gene_counts.txt",
#                     row.names = 1)

cat("Raw counts matrix dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")

# -----------------------------------------------------------------
# STEP 1: TMM normalization with edgeR
# -----------------------------------------------------------------
# Create DGEList object
dge <- DGEList(counts = counts)

# Calculate TMM normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# Compute TMM-normalized counts per million (CPM)
# (prior.count = 0 to match GTEx; log = FALSE for raw CPM)
tmm_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)

cat("Step 1 complete: TMM normalization applied.\n")

# -----------------------------------------------------------------
# STEP 2: Autosomal gene filtering
#
# Requires BOTH conditions met in ≥20% of samples:
#   a) ≥ 0.1 TPM  (using TMM-normalized CPM as TPM proxy, which is
#      the GTEx approach when gene lengths are not available)
#   b) ≥ 6 raw reads (unnormalized counts)
#
# NOTE: If you have gene lengths available, replace tmm_cpm with
#       true TPM values computed as:
#         rpk  <- sweep(counts, 1, gene_lengths_kb, "/")
#         tpm  <- sweep(rpk,    2, colSums(rpk) / 1e6, "/")
#       and use tpm in place of tmm_cpm for the TPM filter below.
# -----------------------------------------------------------------

# --- 2a. Identify autosomal genes ---
# Assumes row names are Ensembl IDs (e.g. ENSG...) or gene symbols.
# Adjust the pattern below to match your annotation.
# Common approach: exclude sex chromosomes (chrX, chrY) and
# mitochondrial (chrM) genes using a gene annotation table.

# Option A (recommended): filter using a BioMart / GTF annotation
#   Example (uncomment and adapt if you have an annotation data frame):
#
# library(biomaRt)
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_info <- getBM(
#   attributes = c("ensembl_gene_id", "chromosome_name"),
#   filters    = "ensembl_gene_id",
#   values     = rownames(counts),
#   mart       = mart
# )
# autosomal_genes <- gene_info$ensembl_gene_id[
#   gene_info$chromosome_name %in% as.character(1:22)
# ]

# Option B (simple fallback): keep all genes if no annotation is
#   available and autosomal filtering has already been applied upstream.
#   Replace the line below with your autosomal gene vector.
autosomal_genes <- fread("/work/jb621/GEMMAmod_iPSCORE/gene_info.txt", data.table = F) %>% filter(chrom %in% paste0("chr", c(1:22)))

# Subset to autosomal genes
tmm_cpm_auto  <- tmm_cpm[rownames(tmm_cpm) %in% autosomal_genes$gene_id, ]
counts_auto   <- counts[rownames(counts) %in% autosomal_genes$gene_id, ]

cat("Autosomal genes retained:", nrow(tmm_cpm_auto), "\n")

# --- 2b. Expression thresholds ---
n_samples     <- ncol(counts_auto)
threshold_pct <- 0.20                          # 20% of samples
min_samples   <- ceiling(n_samples * threshold_pct)

# Filter criterion 1: ≥ 0.1 TPM in ≥ 20% of samples
tpm_filter  <- rowSums(tmm_cpm_auto >= 0.1) >= min_samples

# Filter criterion 2: ≥ 6 raw reads in ≥ 20% of samples
read_filter <- rowSums(counts_auto  >= 6)   >= min_samples

# Both criteria must be satisfied
keep        <- tpm_filter & read_filter

tmm_cpm_filtered <- tmm_cpm_auto[keep, ]
cat("Genes passing both filters:", nrow(tmm_cpm_filtered), "\n")
cat("Step 2 complete: Autosomal gene filtering applied.\n")

# -----------------------------------------------------------------
# STEP 3: Inverse normal transformation (INT) across samples
#         Applied per gene (row-wise across columns/samples)
# -----------------------------------------------------------------
# INT formula (same as GTEx / rank-based):
#   x_int = qnorm( (rank(x) - 0.5) / N )
# where N = number of non-NA observations

inverse_normal_transform <- function(x) {
  # Ties method "average" matches GTEx behavior
  qnorm((rank(x, na.last = "keep", ties.method = "average") - 0.5) /
          sum(!is.na(x)))
}

# Apply row-wise (each gene transformed across samples)
expr_int <- t(apply(tmm_cpm_filtered, 1, inverse_normal_transform))
colnames(expr_int) <- colnames(tmm_cpm_filtered)

cat("Step 3 complete: Inverse normal transformation applied.\n")
cat("Final expression matrix:", nrow(expr_int), "genes x",
    ncol(expr_int), "samples\n")

# -----------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------
# Normalized + filtered + INT-transformed matrix (for eQTL input)
write.table(expr_int,
            file      = "~/ppc_expr_int_normalized.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = NA)   # col.names = NA writes row names correctly

# TMM CPM (pre-INT, post-filter) — useful for QC / visualization
write.table(tmm_cpm_filtered,
            file      = "~/ppc_tmm_cpm_filtered.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = NA)

cat("\nOutputs written:\n")
cat("  ppc_expr_int_normalized.txt  — eQTL-ready INT expression matrix\n")
cat("  ppc_tmm_cpm_filtered.txt     — TMM CPM (pre-INT, for QC)\n")

## =========================================
## VOOM
## =========================================
library(limma)

dge <- DGEList(counts = counts)  # already autosomal-filtered
dge <- calcNormFactors(dge, method = "TMM")

# min.count=10 and min.prop=0.7 are defaults — tune to match
# your desired stringency. To approximate the original criteria:
keep <- filterByExpr(dge, min.count = 6, min.prop = 0.20)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

# If you have a design matrix (e.g. for covariates):
v <- voom(dge_filtered, design, plot = TRUE)


## ============================================================
## ATAC-seq TMM Normalization & Filtering Pipeline
## ============================================================
## Steps:
##   1. TMM normalization (edgeR)
##   2. Remove sex chromosome peaks
##   3. Filter low-accessibility peaks
##      (TMM CPM >= 1 in >= 20% of samples)
## ============================================================

# -----------------------------------------------------------------
# 0. Load data
# -----------------------------------------------------------------
atac_counts <- read.delim("~/ppc_atac_counts.txt",
                          row.names = 1)

cat("Raw ATAC counts matrix dimensions:",
    nrow(atac_counts), "peaks x", ncol(atac_counts), "samples\n")

# Need to map atac_peaks onto chromosome locations so I can filter out the sex chr
atac_peaks <- read.delim("~/ppc_atac_peak.bed")
sex_chroms <- c("chrX", "chrY", "chrM")
atac_peaks_auto <- atac_peaks[!atac_peaks$Chromosome %in% sex_chroms, ]

# -----------------------------------------------------------------
# STEP 1: TMM normalization with edgeR
# -----------------------------------------------------------------
dge <- DGEList(counts = atac_counts)
dge <- calcNormFactors(dge, method = "TMM")

# TMM-normalized CPM (not log-transformed)
tmm_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)

cat("Step 1 complete: TMM normalization applied.\n")

# -----------------------------------------------------------------
# STEP 2: Remove sex chromosome peaks
#
# Assumes row names are in one of these common formats:
#   "chrX:1234-5678"  (UCSC-style, colon-separated)
#   "chrX_1234_5678"  (underscore-separated)
#   "chrX-1234-5678"  (dash-separated)
#
# The regex below extracts the chromosome field (first token)
# and excludes chrX, chrY, and chrM peaks.
# Adjust the pattern if your peak IDs use a different delimiter.
# -----------------------------------------------------------------

peak_names <- atac_peaks_auto$Element_ID

tmm_cpm_auto    <- tmm_cpm[peak_names, ]
atac_counts_auto <- atac_counts[peak_names, ]

cat("Peaks after removing sex chromosomes:", nrow(tmm_cpm_auto), "\n")
cat("Step 2 complete: Sex chromosome peaks removed.\n")

# -----------------------------------------------------------------
# STEP 3: Filter low-accessibility peaks
#
# Keep peaks with TMM CPM >= 1 in >= 20% of samples
# -----------------------------------------------------------------

n_samples    <- ncol(tmm_cpm_auto)
min_samples  <- ceiling(n_samples * 0.20)   # 20% threshold

accessibility_filter <- rowSums(tmm_cpm_auto >= 1) >= min_samples

tmm_cpm_filtered <- tmm_cpm_auto[accessibility_filter, ]

cat("Peaks passing accessibility filter (TMM >= 1 in >= 20% samples):",
    nrow(tmm_cpm_filtered), "\n")
cat("Step 3 complete: Low-accessibility peaks removed.\n")

cat("\nFinal ATAC matrix:", nrow(tmm_cpm_filtered), "peaks x",
    ncol(tmm_cpm_filtered), "samples\n")

# -----------------------------------------------------------------
# STEP 4: Diagnostics — voom trend, zero-inflation, sample outliers,
#          and GEMMAmod lambda deflation check
# -----------------------------------------------------------------
library(limma)
library(ggplot2)
library(patchwork)

# ── 4a. voom mean-variance trend (before & after filtering) ──────
# Before: all autosomal peaks
dge_before <- DGEList(counts = atac_counts_auto)
dge_before <- calcNormFactors(dge_before, method = "TMM")

atac_counts_filtered <- atac_counts_auto[rownames(tmm_cpm_filtered),]

# After: filtered peaks — rebuild DGEList from raw counts so voom
# recomputes library sizes correctly on the filtered set
dge_after <- DGEList(counts = atac_counts_filtered)
dge_after <- calcNormFactors(dge_after, method = "TMM")

par(mfrow = c(1, 2))
v_before <- voom(dge_before, plot = TRUE)
v_after  <- voom(dge_after,  plot = TRUE)
cat("Step 4a complete: voom mean-variance plots generated.\n")

# ── 4b. Zero-inflation in low-count peaks ────────────────────────
log_mean_cpm  <- rowMeans(log2(tmm_cpm_filtered + 0.5))
low_peak_mask <- log_mean_cpm < 2
low_peaks     <- rownames(tmm_cpm_filtered)[low_peak_mask]

cat("\nPeaks with log2(mean TMM-CPM + 0.5) < 2:", sum(low_peak_mask), "\n")

if (length(low_peaks) > 0) {
  zero_frac <- rowMeans(atac_counts_filtered[low_peaks, ] == 0)
  
  par(mfrow = c(1, 2))
  hist(zero_frac, breaks = 50, col = "steelblue",
       main = "Zero fraction in low-count peaks",
       xlab = "Proportion of samples with zero counts")
  
  # Compare strict filter (30%) vs current (20%)
  min_samples_strict   <- ceiling(n_samples * 0.30)
  keep_strict          <- rowSums(tmm_cpm_auto >= 1) >= min_samples_strict
  tmm_cpm_strict       <- tmm_cpm_auto[keep_strict, ]
  atac_counts_strict   <- atac_counts_auto[keep_strict, ]
  
  cat("Peaks removed by stricter 30% filter:",
      nrow(tmm_cpm_filtered) - nrow(tmm_cpm_strict), "\n")
  
  dge_strict  <- DGEList(counts = atac_counts_strict)
  dge_strict  <- calcNormFactors(dge_strict, method = "TMM")
  v_strict    <- voom(dge_strict, plot = TRUE)
} else {
  cat("No low-count peaks remaining — zero inflation not a concern.\n")
}
cat("Step 4b complete: Zero-inflation assessed.\n")

# ── 4c. Sample outlier detection ─────────────────────────────────
log_cpm_mat <- log2(tmm_cpm_filtered + 0.5)

# Library sizes
par(mfrow = c(1, 2))
lib_sizes <- colSums(atac_counts_filtered) / 1e6
barplot(lib_sizes,
        main  = "Library sizes (millions)",
        ylab  = "Total counts (M)",
        col   = ifelse(lib_sizes < median(lib_sizes) * 0.5, "red", "steelblue"),
        las   = 2,
        cex.names = 0.6)
abline(h = median(lib_sizes), col = "red", lty = 2, lwd = 2)
legend("topright", legend = "Median", col = "red", lty = 2, bty = "n")

# PCA
pca     <- prcomp(t(log_cpm_mat), scale. = TRUE)
pca_var <- round(summary(pca)$importance[2, 1:2] * 100, 1)
pca_df  <- data.frame(
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2],
  Sample = colnames(log_cpm_mat)
)

# Flag potential outliers as >3 SD from mean on PC1 or PC2
pca_df$outlier <- abs(scale(pca_df$PC1)) > 3 | abs(scale(pca_df$PC2)) > 3
n_outliers      <- sum(pca_df$outlier)
cat("\nPotential PCA outliers (>3 SD on PC1 or PC2):", n_outliers, "\n")
if (n_outliers > 0) cat("Outlier samples:", pca_df$Sample[pca_df$outlier], "\n")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = outlier, label = Sample)) +
  geom_point(size = 2.5) +
  ggrepel::geom_text_repel(data = pca_df[pca_df$outlier, ],
                           size = 3, max.overlaps = 20) +
  scale_colour_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
  labs(title  = "PCA of log-CPM ATAC-seq counts",
       x      = paste0("PC1 (", pca_var[1], "%)"),
       y      = paste0("PC2 (", pca_var[2], "%)"),
       colour = "Outlier (>3 SD)") +
  theme_minimal(base_size = 12)
print(p_pca)
cat("Step 4c complete: Sample outlier detection done.\n")

# ── 4d. Genomic inflation (lambda) for GEMMAmod caQTL results ────
# Requires gemmamod_caqtl_results with a p_wald column.
# Swap in your actual caQTL results object if named differently.

compute_lambda <- function(p_vals) {
  p_vals  <- p_vals[!is.na(p_vals) & p_vals > 0 & p_vals <= 1]
  chisq   <- qchisq(1 - p_vals, df = 1)
  round(median(chisq) / qchisq(0.5, df = 1), 3)
}

if (exists("gemmamod_caqtl_results") && exists("gemma_caqtl_results")) {
  
  lambda_gemma    <- compute_lambda(gemma_caqtl_results$p_wald)
  lambda_gemmamod <- compute_lambda(gemmamod_caqtl_results$p_wald)
  
  cat("\n--- Genomic inflation (lambda) ---\n")
  cat("Lambda – GEMMA    :", lambda_gemma,    "\n")
  cat("Lambda – GEMMAmod :", lambda_gemmamod, "\n")
  if (lambda_gemmamod < 0.9)
    cat("WARNING: GEMMAmod lambda < 0.9 — likely over-correction by kinship matrix.\n")
  if (lambda_gemmamod > 1.1)
    cat("WARNING: GEMMAmod lambda > 1.1 — residual population structure present.\n")
  
  # QQ plots side by side
  par(mfrow = c(1, 2))
  for (nm in c("GEMMA", "GEMMAmod")) {
    p_vals <- if (nm == "GEMMA") gemma_caqtl_results$p_wald else gemmamod_caqtl_results$p_wald
    p_vals <- sort(p_vals[!is.na(p_vals)])
    n      <- length(p_vals)
    expected <- -log10(seq_len(n) / (n + 1))
    observed <- -log10(p_vals)
    lam      <- compute_lambda(p_vals)
    plot(expected, observed, pch = 19, cex = 0.3,
         col  = ifelse(nm == "GEMMA", "#2C7BB6", "#D7191C"),
         main = paste0("QQ plot – ", nm, "\nλ = ", lam),
         xlab = "Expected –log10(p)",
         ylab = "Observed –log10(p)")
    abline(0, 1, col = "black", lwd = 1.5)
  }
  
  # p-value histograms
  par(mfrow = c(1, 2))
  hist(gemma_caqtl_results$p_wald, breaks = 50, col = "#2C7BB6",
       main = "p-value distribution – GEMMA",
       xlab = "p_wald")
  hist(gemmamod_caqtl_results$p_wald, breaks = 50, col = "#D7191C",
       main = "p-value distribution – GEMMAmod",
       xlab = "p_wald")
  
} else {
  cat("\nSkipping Step 4d: caQTL result objects not found.\n")
  cat("Re-run after loading gemma_caqtl_results and gemmamod_caqtl_results.\n")
}
cat("Step 4d complete: Lambda and p-value distribution checks done.\n")

# ── 4e. voom precision weights check ─────────────────────────────
cat("\n--- voom precision weights summary ---\n")
cat("E matrix dimensions    :", nrow(v_after$E), "x", ncol(v_after$E), "\n")
cat("Weights matrix dimensions:", nrow(v_after$weights), "x", ncol(v_after$weights), "\n")
cat("Weights range          :", round(range(v_after$weights), 3), "\n")
cat("Weights SD             :", round(sd(as.vector(v_after$weights)), 3),
    "(should be > 0; if ~0, weights are constant and not being used)\n")

par(mfrow = c(1, 1))
hist(as.vector(v_after$weights), breaks = 100, col = "steelblue",
     main = "voom precision weights distribution (filtered peaks)",
     xlab = "Weight")
cat("Step 4e complete: voom weights verified.\n")

dev.off()
cat("\nAll diagnostics written to atac_gemmamod_diagnostics.pdf\n")

# Confirm weights are non-trivial
cat("Weights range:", round(range(v_after$weights), 3), "\n")
cat("Weights CV (should be >> 0):", 
    round(sd(v_after$weights) / mean(v_after$weights), 3), "\n")

# -----------------------------------------------------------------
# Save output
# -----------------------------------------------------------------
write.table(tmm_cpm_filtered,
            file      = "~/ppc_atac_tmm_filtered.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = NA)   # col.names = NA writes row names correctly

cat("\nOutput written:\n")
cat("  ppc_atac_tmm_filtered.txt — TMM-normalized, filtered ATAC peak matrix\n")

## =====================
## VOOM ATAC-SEQ
## =====================
dge <- DGEList(counts = atac_counts_auto)  # already sex-chrom filtered
dge <- calcNormFactors(dge, method = "TMM")

# filterByExpr uses a min.count=10 and min.prop=0.7 by default,
# but you can tune these to approximate the original 20% criterion
keep <- filterByExpr(dge, min.count = 1, min.prop = 0.20)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

# voom computes log-CPM and estimates mean-variance weights internally
v <- voom(dge_filtered, design, plot = TRUE)

# Kinship matrix
# To account for genetic relatedness between samples, we performed LD pruning on 
# the 5,536,303 variants using plink83 1.90b6.21 (-indep-pairwise 50 5 0.2). 
# We then used the 323,697 LD-pruned variants to construct a kinship matrix for 
# the 273 iPSCORE individuals using plink83 1.90b6.21 (–make-rel square).

# =============================================================================
# PEER Factor Calculation for PPC eQTL and caQTL
# Replicating: ~25% of samples as max PEER factors, top 2000 most variable features
# PPC eQTL: 22 PEER factors, PPC caQTL: 20 PEER factors
# =============================================================================

library(peer)
library(data.table)

# =============================================================================
# Parameters
# =============================================================================

N_TOP_FEATURES    <- 2000     # top most variable genes/peaks
PEER_PCT_SAMPLES  <- 0.25     # max PEER factors = 25% of samples
set.seed(5366)

# =============================================================================
# Helper: Run PEER on a normalized expression/accessibility matrix
# Input:  expr_mat = samples x features matrix (rows = samples, cols = genes/peaks)
# Output: dataframe of PEER factors (samples x n_factors), colnames = peer1, peer2, ...
# =============================================================================

calculate_peer_factors <- function(expr_mat, n_peer_max) {
  
  message(paste("Input matrix dimensions (samples x features):", 
                nrow(expr_mat), "x", ncol(expr_mat)))
  message(paste("Calculating", n_peer_max, "PEER factors"))
  
  # Select top 2000 most variable features
  feature_vars  <- apply(expr_mat, 2, var)
  top_features  <- names(sort(feature_vars, decreasing = TRUE))[1:min(N_TOP_FEATURES, ncol(expr_mat))]
  expr_subset   <- expr_mat[, top_features]
  
  message(paste("Selected top", ncol(expr_subset), "most variable features"))
  
  # Initialize PEER model
  model <- PEER()
  
  # Set expression matrix (samples x features)
  PEER_setPhenoMean(model, as.matrix(expr_subset))
  
  # Set max number of PEER factors (~25% of samples)
  PEER_setNk(model, n_peer_max)
  
  # Include mean as covariate (recommended)
  PEER_setAdd_mean(model, TRUE)
  
  # Run inference
  PEER_update(model)
  
  # Extract factor matrix (samples x factors)
  factors           <- as.data.frame(PEER_getX(model))
  # Drop the mean factor column (last column when add_mean=TRUE)
  factors           <- factors[, 1:n_peer_max]
  colnames(factors) <- paste0("peer", 1:n_peer_max)
  rownames(factors) <- rownames(expr_mat)
  
  message(paste("PEER factors calculated:", ncol(factors)))
  message(paste("Samples:", nrow(factors)))
  
  return(factors)
}

# =============================================================================
# Load PPC RNA-seq normalized expression matrix
# Expected format: genes x samples (will be transposed to samples x genes)
# =============================================================================

message("Loading PPC RNA-seq normalized expression matrix...")
ppc_rna_mat <- expr_int

# Transpose to samples x genes
ppc_rna_mat <- as.matrix(t(ppc_rna_mat))

n_samples_rna  <- nrow(ppc_rna_mat)
n_peer_max_rna <- floor(n_samples_rna * PEER_PCT_SAMPLES)

message(paste("PPC RNA-seq samples:", n_samples_rna))
message(paste("Max PEER factors (25% of samples):", n_peer_max_rna))

# Calculate PEER factors for eQTL
ppc_eqtl_peer <- calculate_peer_factors(ppc_rna_mat, n_peer_max_rna)

# Save all PEER factors (you will later select optimal n=22 for PPC eQTL)
fwrite(ppc_eqtl_peer, 
       file      = "~/data/PPC_eQTL_peer_factors.txt", 
       sep       = "\t", 
       row.names = TRUE, 
       col.names = TRUE)
message("Saved: PPC_eQTL_peer_factors.txt")

# =============================================================================
# Load PPC ATAC-seq normalized accessibility matrix
# Expected format: peaks x samples (will be transposed to samples x peaks)
# =============================================================================

# message("Loading PPC ATAC-seq normalized accessibility matrix...")
# ppc_atac_mat <- add_rownames(fread("/work/jb621/GEMMAmod_iPSCORE/ppc_atac_peak.bed", data.table = FALSE))
# 
# # Drop genomic coordinate columns if present (chr, start, end)
# # Adjust column indices based on your actual file format
# coord_cols   <- c("chrom", "start", "end")
# coord_cols   <- coord_cols[coord_cols %in% colnames(ppc_atac_mat)]
# if (length(coord_cols) > 0) {
#   ppc_atac_mat <- ppc_atac_mat[, !colnames(ppc_atac_mat) %in% coord_cols]
# }

# Transpose to samples x peaks
ppc_atac_mat <- as.matrix(t(tmm_cpm_filtered))

n_samples_atac  <- nrow(ppc_atac_mat)
n_peer_max_atac <- floor(n_samples_atac * PEER_PCT_SAMPLES)

message(paste("PPC ATAC-seq samples:", n_samples_atac))
message(paste("Max PEER factors (25% of samples):", n_peer_max_atac))

# Calculate PEER factors for caQTL
ppc_caqtl_peer <- calculate_peer_factors(ppc_atac_mat, n_peer_max_atac)

# Save all PEER factors (you will later select optimal n=20 for PPC caQTL)
fwrite(ppc_caqtl_peer, 
       file      = "~/data/PPC_caQTL_peer_factors.txt", 
       sep       = "\t", 
       row.names = TRUE, 
       col.names = TRUE)
message("Saved: PPC_caQTL_peer_factors.txt")

# =============================================================================
# Extract final covariates using the published optimal PEER factor numbers
# PPC eQTL:  22 PEER factors
# PPC caQTL: 20 PEER factors
# =============================================================================

ppc_eqtl_peer_final  <- ppc_eqtl_peer[,  1:22]
ppc_caqtl_peer_final <- ppc_caqtl_peer[, 1:20]

fwrite(ppc_eqtl_peer_final,  
       file      = "/work/jb621/GEMMAmod_iPSCORE/PPC_eQTL_peer_factors_final_n22.txt",  
       sep       = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

fwrite(ppc_caqtl_peer_final, 
       file      = "/work/jb621/GEMMAmod_iPSCORE/PPC_caQTL_peer_factors_final_n20.txt", 
       sep       = "\t", 
       row.names = TRUE, 
       col.names = TRUE)

message("Done. Final PEER factor files saved.")
message("PPC eQTL:  PPC_eQTL_peer_factors_final_n22.txt  (22 factors)")
message("PPC caQTL: PPC_caQTL_peer_factors_final_n20.txt (20 factors)")

# voom version on filtered but not normalized data

# Process covariates
# scale all covariates
scale_metadata = function(metadata_sample) {
  to_scale                   = colnames(metadata_sample)[ !colnames(metadata_sample) %in% c("phenotype_id", "genotype_id", "subject_id")]
  metadata_sample[,to_scale] = scale(metadata_sample[,to_scale])
  
  return(metadata_sample)
}

process_covariates = function(config, peerdata) {
  metadata_sample      = fread(config$metadata_sample , sep = "\t", header = TRUE, data.table = FALSE)
  metadata_subject     = fread(config$metadata_subject, sep = "\t", header = TRUE, data.table = FALSE)
  metadata             = metadata_sample[,c("phenotype_id", "genotype_id", "subject_id")]
  covariates           = merge(scale_metadata(metadata_sample), peerdata, by.x = "phenotype_id", by.y = "row.names")
  covariates           = merge(covariates                     , metadata_subject, by = c("phenotype_id", "genotype_id", "subject_id"))
  rownames(covariates) = covariates$phenotype_id
  rownames(metadata  ) = metadata  $phenotype_id
  
  message(paste("Number of rows in covariates:", nrow(covariates)))
  message(paste("Number of rows in metadata:", nrow(metadata)))
  
  covariates[,c("subject_id", "genotype_id")] = NULL
  
  covariates = covariates[ metadata$phenotype_id,]    
  
  fwrite(covariates, file = paste(config$out_folder, "step_1", "covariates.txt", sep = "/"), sep = "\t", row.names = FALSE, col.names = TRUE)
  fwrite(metadata  , file = paste(config$out_folder, "step_1", "metadata.txt"  , sep = "/"), sep = "\t", row.names = FALSE, col.names = TRUE)
  
  message(paste("Saved:", paste(config$out_folder, "step_1", "covariates.txt", sep = "/")))
  message(paste("Saved:", paste(config$out_folder, "step_1", "metadata.txt"  , sep = "/")))
  
  #return(sort(unique(metadata$genotype_id)))
}

peer_file = paste(config$out_folder, "step_1", "phenotype", paste("peer", "factors"  , "txt", sep = "."), sep = "/")
peerdata  = add_rownames(fread(peer_file, data.table = F))
process_covariates(config, peerdata)


# OPTIONAL
# Aligning caQTL, eQTL, and GWA Intervals 
GWAS_QTL_Colocalization_Summaries <- read_tsv(file = "~/GWAS_QTL_Colocalization_Summaries.txt")
GWAS_QTL_Colocalization_Summaries_PPC <- GWAS_QTL_Colocalization_Summaries %>% 
  dplyr::filter(Tissue == "PPC")

ppc_atac_peak_244298_summary <- GWAS_QTL_Colocalization_Summaries_PPC %>% 
  dplyr::filter(Element_ID == "ppc_atac_peak_244298")

ppc_atac_peak_244305_summary <- GWAS_QTL_Colocalization_Summaries_PPC %>% 
  dplyr::filter(Element_ID == "ppc_atac_peak_244305")

# Cluster is PPC 122
Cluster_PPC_122 <- GWAS_QTL_Colocalization_Summaries_PPC %>% 
  dplyr::filter(Cluster_ID == "PPC_122")

# ENSG00000153814.13 is the Ensembl name for JAZF1 

PPC_ATAC_bed_intervals <- read_tsv("/work/jb621/GEMMAmod_iPSCORE/ppc_atac_peak.bed")
PPC_ATAC_bed_tested_intervals <- PPC_ATAC_bed_intervals %>% 
  dplyr::filter(Expressed == "TRUE")

PPC_ATAC_bed_enriched_intervals <- PPC_ATAC_bed_tested_intervals %>% 
  dplyr::filter(Element_ID %in% c("ppc_atac_peak_244298", "ppc_atac_peak_244305"))

# various QTL results
PPC_caQTLs_sumStats <- fread(file = "~/PPC_discovery_caqtls_stats.txt")
PPC_caQTLs_sumStats_JAZF1 <- PPC_caQTLs_sumStats %>% 
  dplyr::filter(ElementID %in% c("ppc_atac_peak_244298", "ppc_atac_peak_244305"))
# Some SNPs have multiple test statistics, and that's because of the "format that they"discovery" 
# process that they were tested in 
# "To identify additional independent QTL associations for a gene or peak 
# (i.e., conditional QTLs), we performed stepwise regression analysis in which 
# we re-performed QTL analysis with the genotype of the lead variant for the QTL 
# as a covariate. We repeated the procedure to discover up to three conditional 
# associations. For each iteration, we performed the two-step procedure described 
# above and considered conditional eQTLs with q-values <0.05 as significant.

PPC_eQTLs_sumStats <- fread(file = "~/PPC_discovery_eqtls_stats.txt")
PPC_eQTLs_sumStats_JAZF1 <- PPC_eQTLs_sumStats %>% 
  dplyr::filter(ElementID %in% "ENSG00000153814.13")


