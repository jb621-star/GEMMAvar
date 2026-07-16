## ============================================================
## 01_write_peak_inputs.R
## Write per-peak RDS files for caQTL SLURM array
## Mirrors 01_write_gene_inputs.R but for ATAC-seq peaks
##
## Key differences from eQTL:
##   - Window: 100kb either side of peak midpoint (not 1Mb)
##   - MHC exclusion: peaks within 100kb of chr6:28,510,120-33,480,577
##   - All peaks tested regardless of Expressed flag
##   - Peak IDs from peaks$Element_ID (e.g. ppc_atac_peak_1)
## ============================================================

library(limma)
library(edgeR)
library(tidyverse)
#library(Matrix)
library(data.table)
library(parallel)
library(GenomicRanges)
setwd("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/")


# ---- Load data ----------------------------------------------------------------
# Adjust these paths/object names to match your environment
#load(file = "edgeR_atac_ordered_v2.Rda")
#atac_expr_ordered <- edgeR_expr_ordered
#dge_aligned <- dge_aligned[rownames(atac_expr_ordered),]
save(genos_mat_filtered, atac_expr_ordered, dge_aligned, atac_cov_ordered, filtered_kinship, peaks,
     file = "allPairs_testing_inputs_caQTL.Rda")

load("allPairs_testing_inputs_caQTL.Rda")   # should contain:

saveRDS(dge_aligned,"caqtl_gene_inputs/dge_atac.rds")
#   genos_mat_filtered   (105 x 5733665)
#   atac_expr_ordered    (peaks x 105, normalized + INT)
#   atac_cov_ordered      (covariate table)
#   filtered_kinship     (105 x 105)
#   peaks                (peak annotation data.table)

# Build covariate matrix
atac_cov_ordered_matrix <- atac_cov_ordered %>%
  select(!c("Subject_UUID", "ATAC_UUID")) %>%
  as.matrix()
W <- cbind(rep(1, nrow(atac_cov_ordered_matrix)), atac_cov_ordered_matrix)
K <- filtered_kinship

# Output directories
dir.create("caqtl_gene_inputs",   showWarnings = FALSE)
dir.create("caqtl_gene_outputs",  showWarnings = FALSE)

#load("atac_dge_filtered.Rda")                # dge_filtered
#dge_lead     <- dge_aligned[, meta_map_ordered$ATAC_UUID]
#saveRDS(dge_lead, "caqtl_gene_inputs/dge_atac.rds")
#atac_expr_ordered<- readRDS(file = "caqtl_gene_inputs/dge_atac.rds")

# Write shared inputs once
saveRDS(W, "caqtl_gene_inputs/W.rds")
saveRDS(K, "caqtl_gene_inputs/K.rds")

peaks <- peaks[peaks$Element_ID %in% rownames(edgeR_expr_ordered),]

# ---- MHC exclusion (100kb window around MHC region) --------------------------
mhc_chr   <- "chr6"
mhc_start <- 28510120 - 100000
mhc_end   <- 33480577 + 100000

mhc_peaks <- peaks[Chromosome == mhc_chr &
                     End   >= mhc_start &
                     Start <= mhc_end,
                   Element_ID]

cat("Peaks excluded in/near MHC region:", length(mhc_peaks), "\n")

peaks_to_test <- peaks[!Element_ID %in% mhc_peaks]
cat("Peaks remaining after MHC exclusion:", nrow(peaks_to_test), "\n")

# ---- Build SNP GRanges -------------------------------------------------------
snp_ids <- colnames(genos_mat_filtered)

snp_gr <- data.frame(snp_id = snp_ids) %>%
  tidyr::separate(snp_id, into = c("chr", "pos", "ref", "alt"),
                  sep = "_", remove = FALSE, extra = "drop") %>%
  mutate(pos = as.integer(pos)) %>%
  with(GRanges(
    seqnames = chr,
    ranges   = IRanges(start = pos, end = pos),
    snp_id   = snp_id
  ))

peak_windows <- peaks_to_test %>%
  mutate(
    win_start = pmax(1, Start),
    win_end   = End
  ) %>%
  with(GRanges(
    seqnames  = Chromosome,
    ranges    = IRanges(start = win_start, end = win_end),
    peak_id   = Element_ID
  ))

# ---- Find overlaps -----------------------------------------------------------
hits <- findOverlaps(peak_windows, snp_gr)

testing_pairs_caqtl <- data.frame(
  Element_ID = peak_windows$peak_id[queryHits(hits)],
  SNP_ID     = snp_gr$snp_id[subjectHits(hits)]
)

cat("Total caQTL testing pairs:", nrow(testing_pairs_caqtl), "\n")
cat("Peaks with >= 1 SNP      :", n_distinct(testing_pairs_caqtl$Element_ID), "\n")
cat("Avg SNPs per peak        :", round(nrow(testing_pairs_caqtl) /
                                          n_distinct(testing_pairs_caqtl$Element_ID)), "\n")

# Save full pairs table for reference
fwrite(testing_pairs_caqtl, "caqtl_gene_inputs/testing_pairs_caqtl.tsv", sep = "\t")
testing_pairs_caqtl <- fread("caqtl_gene_inputs/testing_pairs_caqtl.tsv")

# ---- Write per-peak RDS files ------------------------------------------------
peak_list <- missing_peaks 
peak_list      <- unique(testing_pairs_caqtl$Element_ID)
peaks_written  <- character(0)
n_missing_snps <- 0

cat("Writing per-peak input files...\n")

for (i in seq_along(peak_list)) {
  
  peak_id <- peak_list[i]
  if (i %% 500 == 0) cat(sprintf("  %d / %d peaks\n", i, length(peak_list)))
  
  snp_ids_peak  <- testing_pairs_caqtl$SNP_ID[testing_pairs_caqtl$Element_ID == peak_id]
  snp_ids_found <- snp_ids_peak[snp_ids_peak %in% colnames(genos_mat_filtered)]
  
  n_missing_snps <- n_missing_snps + (length(snp_ids_peak) - length(snp_ids_found))
  if (length(snp_ids_found) == 0) next
  
  # Genotype matrix for this peak: samples x SNPs
  geno_sub <- genos_mat_filtered[, snp_ids_found, drop = FALSE]
  
  # Mean-impute missing genotypes
  for (j in seq_len(ncol(geno_sub))) {
    nas <- is.na(geno_sub[, j])
    if (any(nas)) geno_sub[nas, j] <- mean(geno_sub[, j], na.rm = TRUE)
  }
  
  # Accessibility vector for this peak
  expr_vec <- as.numeric(atac_expr_ordered[peak_id, ])
  
  saveRDS(
    list(peak_id  = peak_id,
         geno     = geno_sub,
         expr     = expr_vec,
         snp_ids  = snp_ids_found),
    file = paste0("caqtl_gene_inputs/", peak_id, ".rds")
  )
  
  peaks_written <- c(peaks_written, peak_id)
}

writeLines(peaks_written, "caqtl_gene_inputs/peak_list.txt")
#writeLines(peak_list, "caqtl_gene_inputs/peak_list.txt")

cat("\nDone.\n")
cat("Peaks written      :", length(peaks_written), "\n")
cat("Missing SNPs total :", n_missing_snps, "\n")
cat("Peak list written to caqtl_gene_inputs/peak_list.txt\n")