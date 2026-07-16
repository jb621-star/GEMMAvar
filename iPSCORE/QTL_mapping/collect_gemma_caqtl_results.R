## ============================================================
## 04_collect_results_caQTL.R
## Collect per-peak output RDS files into a single results table
## Run after all SLURM jobs complete
## ============================================================

library(data.table)

peak_list <- readLines("caqtl_gene_inputs/peak_list.txt")
cat("Peaks expected:", length(peak_list), "\n")

# ---- Check completion status ----
output_files <- paste0("caqtl_gemmamod_outputs/", peak_list, "NOSHRINKAGE.rds")
completed    <- file.exists(output_files)
error_files  <- list.files("caqtl_gemmamod_outputs", pattern = "^ERROR_", full.names = TRUE)

cat("Completed :", sum(completed), "\n")
cat("Missing   :", sum(!completed), "\n")
cat("Errors    :", length(error_files), "\n")

if (sum(!completed) > 0) {
  missing_peaks <- peak_list[!completed]
  cat("\nFirst 20 incomplete peaks:\n")
  print(head(missing_peaks, 20))
  writeLines(missing_peaks, "caqtl_gene_inputs/peak_list_rerun.txt")
  cat("Full list written to caqtl_gene_inputs/peak_list_rerun.txt\n")
  cat("Resubmit with:\n")
  cat("  sbatch --array=1-", length(missing_peaks), "%1500 03_submit_caqtl_array.sh\n",
      sep = "")
  cat("  (after replacing peak_list.txt with peak_list_rerun.txt)\n")
}

if (length(error_files) > 0) {
  cat("\nError details:\n")
  for (f in error_files) {
    err <- readRDS(f)
    cat(sprintf("  %s: %s\n", err$peak_id, err$error))
  }
}

# ---- Collect completed results ----
cat("\nCollecting results...\n")

results_list <- vector("list", sum(completed))
j <- 0

for (i in which(completed)) {
  peak_id <- peak_list[i]
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

# Standardise column names
setnames(results,
         old = c("rs"),
         new = c("SNP_ID"),
         skip_absent = TRUE)

col_order <- c("SNP_ID", "Element_ID", "beta", "se_beta",
               "tau", "lambda", "F_wald", "p_wald")
col_order <- col_order[col_order %in% colnames(results)]
setcolorder(results, col_order)

cat(sprintf("Total pairs collected : %d\n", nrow(results)))
cat(sprintf("Peaks collected       : %d\n", uniqueN(results$Element_ID)))

fwrite(results, "gemmamod_caqtl_results_NOSHRINKAGE_allPairs.tsv", sep = "\t")
save(results,   file = "gemmamod_caqtl_results_NOSHRINKAGE_allPairs.Rda")
cat("Results written to gemmamod_caqtl_results_NOSHRINKAGE_allPairs.tsv and .Rda\n")