## ============================================================
## 04_collect_results.R
## Collect per-gene output RDS files into a single results table
## Run after all SLURM jobs complete
## ============================================================

library(data.table)

setwd("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping")

gene_list <- readLines("gemma_gene_inputs/gene_list.txt")
cat("Genes expected:", length(gene_list), "\n")

# ---- Check job completion status ----
output_files <- paste0("gemmamod_gene_outputs/", gene_list, ".rds")
completed    <- file.exists(output_files)
error_files  <- list.files("gemmamod_gene_outputs", pattern = "^ERROR_", full.names = TRUE)

cat("Completed :", sum(completed), "\n")
cat("Missing   :", sum(!completed), "\n")
cat("Errors    :", length(error_files), "\n")

if (sum(!completed) > 0) {
  cat("\nGenes not yet completed:\n")
  print(gene_list[!completed])
  cat("\nRe-submit missing jobs before collecting, or proceed with partial results.\n")
}

if (length(error_files) > 0) {
  cat("\nError details:\n")
  for (f in error_files) {
    err <- readRDS(f)
    cat(sprintf("  %s: %s\n", err$gene_id, err$error))
  }
}

# In R: write the incomplete genes to a rerun list
incomplete_genes <- gene_list[!completed]
writeLines(incomplete_genes, "gemma_gene_inputs/gene_list_rerun.txt")
cat("Wrote", length(incomplete_genes), "genes to gene_list_rerun.txt\n")

# ---- Collect completed results ----
cat("\nCollecting results...\n")

results_list <- vector("list", sum(completed))
j <- 0

for (i in which(completed)) {
  gene_id <- gene_list[i]
  res <- tryCatch(
    readRDS(output_files[i]),
    error = function(e) {
      cat("Failed to read:", output_files[i], "\n"); NULL
    }
  )
  if (is.null(res)) next
  j <- j + 1
  results_list[[j]] <- as.data.table(res)
}

results <- rbindlist(results_list[seq_len(j)], fill = TRUE)

# Standardise column names to match your existing def_results format
setnames(results,
         old = c("rs"),          # pygemma's default SNP column name
         new = c("SNP_ID"),
         skip_absent = TRUE)

# Reorder columns
col_order <- c("SNP_ID", "Element_ID", "beta", "se_beta",
               "tau", "lambda", "F_wald", "p_wald")
col_order <- col_order[col_order %in% colnames(results)]
setcolorder(results, col_order)

cat(sprintf("\nTotal pairs collected: %d\n", nrow(results)))
cat(sprintf("Genes collected      : %d\n", uniqueN(results$Element_ID)))

# Save
fwrite(results, "gemmamod_results_allPairs.tsv", sep = "\t")
save(results,   file = "gemmamod_results_allPairs.Rda")
cat("Results written to gemmamod_results_allPairs.tsv and .Rda\n")
