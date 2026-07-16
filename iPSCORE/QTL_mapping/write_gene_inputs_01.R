## ============================================================
## 01_write_gene_inputs.R
## Write per-gene RDS files containing SNPs, phenotype, W, K
## Run once interactively before submitting SLURM array
## ============================================================
.libPaths(c("/home/jb621/R/x86_64-pc-linux-gnu-library/4.2","/usr/local/lib/R/site-library","/usr/local/lib/R/library",.libPaths()) ) 
library(limma)
library(edgeR)
library(tidyverse)
#library(Matrix)
library(data.table)
library(parallel)

setwd("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/")

load("allPairs_testing_inputs_050226.Rda")

# Build covariate matrix (same as your existing code)
rna_cov_ordered_matrix <- rna_cov_ordered %>%
  select(!c("Subject_UUID", "RNA_UUID")) %>%
  as.matrix()
W <- cbind(rep(1, nrow(rna_cov_ordered_matrix)), rna_cov_ordered_matrix)
K <- filtered_kinship

# Output directories
dir.create("gemma_gene_inputs",  showWarnings = FALSE)
dir.create("gemma_gene_outputs", showWarnings = FALSE)

# Get unique genes from testing_pairs (already 1Mb-windowed, MHC-excluded)
gene_list <- unique(testing_pairs$Element_ID)
cat("Total genes to process:", length(gene_list), "\n")

#write_tsv(data.frame(gene_list), file = "gemma_gene_inputs/gene_list.txt", col_names = F)

# Write W and K once (shared across all jobs)
saveRDS(W, "gemma_gene_inputs/W.rds")
saveRDS(K, "gemma_gene_inputs/K.rds")
saveRDS(dge_lead, "gemma_gene_inputs/dge_lead.rds")

# Write per-gene RDS: SNP genotypes + expression vector
n_missing_snps <- 0
genes_written  <- character(0)

#for (i in 14000:length(gene_list)) {
for (i in seq_along(gene_list)) {
  
  gene_id <- gene_list[i]
  if (i %% 500 == 0) cat(sprintf("  %d / %d genes\n", i, length(gene_list)))
  
  # SNPs to test for this gene
  snp_ids <- testing_pairs$SNP_ID[testing_pairs$Element_ID == gene_id]
  
  # Filter to SNPs present in genotype matrix
  snp_ids_found <- snp_ids[snp_ids %in% colnames(genos_mat_filtered)]
  n_missing_snps <- n_missing_snps + (length(snp_ids) - length(snp_ids_found))
  
  if (length(snp_ids_found) == 0) next
  
  # Genotype matrix for this gene: samples x SNPs
  geno_sub <- genos_mat_filtered[, snp_ids_found, drop = FALSE]
  
  # Mean-impute any missing genotypes
  for (j in seq_len(ncol(geno_sub))) {
    nas <- is.na(geno_sub[, j])
    if (any(nas)) geno_sub[nas, j] <- mean(geno_sub[, j], na.rm = TRUE)
  }
  
  # Expression vector for this gene
  expr_vec <- as.numeric(edgeR_expr_ordered[gene_id, ])
  
  saveRDS(
    list(gene_id = gene_id,
         geno    = geno_sub,
         expr    = expr_vec,
         snp_ids = snp_ids_found),
    file = paste0("gemma_gene_inputs/", gene_id, ".rds")
  )
  
  genes_written <- c(genes_written, gene_id)
}

# Write gene list for SLURM array indexing
writeLines(genes_written, "gemma_gene_inputs/gene_list.txt")
gene_list <- readLines("gemma_gene_inputs/gene_list.txt")

cat("\nDone.\n")
cat("Genes written      :", length(genes_written), "\n")
cat("Missing SNPs total :", n_missing_snps, "\n")
cat("Gene list written to gemma_gene_inputs/gene_list.txt\n")