.libPaths(c("/home/jb621/R/x86_64-pc-linux-gnu-library/4.4","/usr/local/lib/R/site-library","/usr/local/lib/R/library",.libPaths()) ) 
library(parallel)
library(coloc)
library(data.table)
library(dplyr)
library(optparse)
library(pbapply)

setwd("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/")

#load(file = "PPC_FDR05_eQTL_SNPs_bothGEMMA.Rda")
#load(file = "PPC_FDR05_caqtl_SNPs_bothGEMMA.Rda")

load(file = "PPC_eQTL_SNPs_bothGEMMA.Rda")
load(file = "PPC_caqtl_SNPs_bothGEMMA.Rda")

# head(gemma_eQTL_results_loc)
# head(gemma_caqtl_results_loc)
# head(gemmamod_eQTL_results_loc)
# head(gemmamod_caqtl_results_loc)
# 
# #save(gemma_eQTL_results_loc, gemmamod_eQTL_results_loc, file = "PPC_eQTL_SNPs_bothGEMMA.Rda")
# #save(gemma_caqtl_results_loc, gemmamod_caqtl_results_loc, file = "PPC_caqtl_SNPs_bothGEMMA.Rda")
# 
# # Check how many SNPs are shared per gene-peak pair for a sample of pairs
# eqtl_gene  <- gemma_sig[Element_ID == gemma_sig$Element_ID[1]]
# caqtl_peak <- gemma_caqtl_sig[Element_ID == gemma_caqtl_sig$Element_ID[1]]
# 
# shared_snps <- intersect(eqtl_gene$SNP, caqtl_peak$SNP)
# cat("Shared SNPs for first gene-peak pair:", length(shared_snps), "\n")
# 
# # Check across a random sample of 100 pairs
# gene_sample <- sample(unique(gemma_sig$Element_ID), 100)
# peak_sample <- sample(unique(gemma_caqtl_sig$Element_ID), 100)
# 
# shared_counts <- sapply(seq_along(gene_sample), function(i) {
#   g <- gemma_sig[Element_ID == gene_sample[i]]
#   p <- gemma_caqtl_sig[Element_ID == peak_sample[i]]
#   length(intersect(g$SNP, p$SNP))
# })
# 
# summary(shared_counts)
# cat("Pairs with >= 50 shared SNPs:", sum(shared_counts >= 50), "\n")
# cat("Pairs with 0 shared SNPs    :", sum(shared_counts == 0), "\n")
# 
# # Pairs sharing at least one significant SNP
# candidate_pairs <- merge(gemma_sig, gemma_caqtl_sig, by = "SNP")
# cat("Candidate gene-peak pairs to test:", nrow(candidate_pairs), "\n")

# ── Load expression matrices once ──────────────────────────────────────────────
load("allPairs_testing_inputs_050226.Rda")   # loads edgeR_expr_ordered
load("allPairs_testing_inputs_caQTL.Rda")    # loads atac_expr_ordered

# ── Build the pairs table ───────────────────────────────────────────────────────
# Option A: if you have a pre-made table with columns gene_id & peak_id
# pairs <- fread("gene_peak_pairs.csv")

# Keep only gene/peak pairs that share at least one SNP
# gemmamod_caqtl_sig <- gemmamod_caqtl_sig %>%
#   dplyr::rename(SNP = SNP_ID)

# gemmamod_sig <- gemmamod_sig %>%
#   dplyr::rename(SNP = SNP_ID)

gemmamod_caqtl_results_loc <- gemmamod_caqtl_results_loc %>%
  dplyr::rename(SNP = SNP_ID)
# # 
# # eqtl_snps  <- unique(gemmamod_sig[, .(gene_id = Element_ID, SNP)])
# # SNPs tested in caQTL (use union of both methods if they differ)
# caqtl_snps <- unique(c(
#   gemma_caqtl_results_loc$SNP,
#   gemmamod_caqtl_results_loc$SNP
# ))
# 
# # Find overlapping SNPs for each eQTL method
# gemma_eQTL_caQTL_SNPs <- gemma_eQTL_results_loc %>%
#   dplyr::filter(SNP %in% caqtl_snps) %>%
#   dplyr::rename(gene_id = Element_ID)
# 
# gemmamod_eQTL_caQTL_SNPs <- gemmamod_eQTL_results_loc %>%
#   dplyr::filter(SNP %in% caqtl_snps) %>%
#   rename(gene_id = Element_ID)
# 
# # Step 1: Get the SNPs tested for each gene (already filtered to caQTL SNPs)
# gene_snps <- gemmamod_eQTL_caQTL_SNPs[, .(gene_id, SNP)]
# 
# # Step 2: Get the SNPs tested for each peak
# peak_snps <- gemmamod_caqtl_results_loc[, .(peak_id = Element_ID, SNP)]
# 
# # Step 3: Join on SNP — this gives every gene-peak pair that shares at least one SNP
# gene_peak_snp <- merge(gene_snps, peak_snps, by = "SNP")
# 
# # Step 3.5 filtering to just cases of >1 SNP per gene-peak pair
# # Count SNPs per gene-peak pair, then filter
# gene_peak_snp_filtered <- gene_peak_snp[
#   gene_peak_snp[, .N, by = .(gene_id, peak_id)][N > 1],
#   on = .(gene_id, peak_id)
# ]

# Check how many pairs remain
# cat("Pairs before filter:", uniqueN(gene_peak_snp[, paste(gene_id, peak_id)]), "\n")
# cat("Pairs after filter:", uniqueN(gene_peak_snp_filtered[, paste(gene_id, peak_id)]), "\n")

# # Step 4: Check what this looks like
# head(gene_peak_snp)
# cat("Unique genes:", uniqueN(gene_peak_snp$gene_id), "\n")
# cat("Unique peaks:", uniqueN(gene_peak_snp$peak_id), "\n")
# cat("Unique pairs:", nrow(unique(gene_peak_snp[, .(gene_id, peak_id)])), "\n")
# cat("Unique SNPs:", uniqueN(gene_peak_snp$SNP), "\n")

# Collapse to unique pairs
# pairs <- unique(gene_peak_snp_filtered[, .(gene_id, peak_id)])

# CHUNK_SIZE <- 107070

# # Split into chunks
# pair_chunks <- split(pairs, ceiling(seq_len(nrow(pairs)) / CHUNK_SIZE))
# n_chunks <- length(pair_chunks)
# cat("Number of chunks:", n_chunks, "\n")
# 
# # Save each chunk as a CSV
# dir.create("pair_chunks", showWarnings = FALSE)
# for (i in seq_along(pair_chunks)) {
#   fwrite(pair_chunks[[i]], 
#          file = sprintf("pair_chunks/pairs_chunk_%04d.csv", i))
# }
# cat("Chunks written to pair_chunks\n")

option_list <- list(
  make_option("--chunk_file", type = "character"),
  make_option("--outfile",    type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Load this chunk's pairs
pairs <- fread(opt$chunk_file)

# gene_id <- pairs$gene_id[1]
# peak_id <- pairs$peak_id[1]

# gemma_caqtl_results_loc <- gemma_caqtl_results_loc %>% 
#   dplyr::rename(Element_ID = peak_id)
# gemmamod_caqtl_results_loc <- gemmamod_caqtl_results_loc %>% 
#   dplyr::rename(Element_ID = peak_id)

# ── Helper: run coloc for one pair ─────────────────────────────────────────────
run_coloc <- function(gene_id, peak_id) {
  
  eqtl_sub  <- gemmamod_eQTL_results_loc[Element_ID == gene_id]
  caqtl_sub <- gemmamod_caqtl_results_loc[Element_ID == peak_id]
  
  merged <- merge(eqtl_sub, caqtl_sub, by = "SNP")
  if (nrow(merged) < 2L) return(NULL)
  
  sdY.x <- sd(t(edgeR_expr_ordered[gene_id, ]))
  sdY.y <- sd(t(atac_expr_ordered[peak_id,  ]))
  
  res <- tryCatch(
      coloc.abf(
        dataset1 = list(
          snp     = merged$SNP,
          beta    = merged$beta.x,
          varbeta = merged$se_beta.x^2,
          type    = "quant",
          sdY     = sdY.x,
          N       = 105
        ),
        dataset2 = list(
          snp     = merged$SNP,
          beta    = merged$beta.y,
          varbeta = merged$se_beta.y^2,
          type    = "quant",
          sdY     = sdY.y,
          N       = 105
        )
      ),
    error = function(e) {
      message("coloc failed for ", gene_id, " / ", peak_id, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(res)) return(NULL)
  
  data.table(
    gene_id = gene_id,
    peak_id = peak_id,
    nsnps   = res$summary["nsnps"],
    PP.H0   = res$summary["PP.H0.abf"],
    PP.H1   = res$summary["PP.H1.abf"],
    PP.H2   = res$summary["PP.H2.abf"],
    PP.H3   = res$summary["PP.H3.abf"],
    PP.H4   = res$summary["PP.H4.abf"]
  )
}

# ── Run over all pairs with progress bar ───────────────────────────────────────
results_list <- pbmapply(
  run_coloc,
  gene_id  = pairs$gene_id,
  peak_id  = pairs$peak_id,
  SIMPLIFY = FALSE
)

results <- rbindlist(Filter(Negate(is.null), results_list))

dir.create("coloc_results", showWarnings = FALSE)
fwrite(results, opt$outfile)
cat("Done:", nrow(results), "results written to", opt$outfile, "\n")

# Read the results and compare
# load(file = "gemma_coloc_results_all_pairs_051626.Rda")
# gemma_results <- results
# 
# load(file = "gemmamod_coloc_results_all_pairs.Rda")
# gemmamod_results <- results
# 
# # GEMMAmod's summary statistics allowed for 27 additional gene/peak pairs for testing
# # Both versions were able to colocalize JAZF1 (ENSG00000153814.13) and ppc_atac_peak_244298
# # GEMMAmod was more certain of this colocalization than GEMMA, but only from 0.99999 to 1
# # Assuming both results are loaded as gemma_coloc and gemmamod_coloc
# 
# # Classify each pair in each result set
# classify_coloc <- function(dt, pp4_thresh = 0.75, pp4_3_thresh = 3) {
#   dt[, coloc := PP.H4 >= pp4_thresh & (PP.H4 / PP.H3) >= pp4_3_thresh]
#   dt
# }
# 
# gemma_coloc    <- classify_coloc(gemma_results)
# gemmamod_coloc <- classify_coloc(gemmamod_results)
# 
# # Merge on gene/peak pair
# comparison <- merge(
#   gemma_coloc[,    .(gene_id, peak_id, PP.H4.gemma    = PP.H4, PP.H3.gemma    = PP.H3, coloc.gemma    = coloc)],
#   gemmamod_coloc[, .(gene_id, peak_id, PP.H4.gemmamod = PP.H4, PP.H3.gemmamod = PP.H3, coloc.gemmamod = coloc)],
#   by = c("gene_id", "peak_id")
# )
# 
# # Four categories
# comparison[, category := fcase(
#   coloc.gemma &  coloc.gemmamod, "shared",
#   coloc.gemma & !coloc.gemmamod, "GEMMA only",
#   !coloc.gemma &  coloc.gemmamod, "GEMMAmod only",
#   !coloc.gemma & !coloc.gemmamod, "neither"
# )]
# 
# table(comparison$category)
# 
# gemma_results_H4H3 <- gemma_results %>% 
#   filter(PP.H4 >= 0.75 & PP.H4 / PP.H3 >= 3)
# # 863 / 2120 gene-peak pairs meet the colocalization criteria 
# 
# gemmamod_results_H4H3 <- gemmamod_results %>% 
#   filter(PP.H4 >= 0.75 & PP.H4 / PP.H3 >= 3)
# # 990 / 2147 gene-peak pairs meet the colocalization criteria 
# 
# gemmamod_results |>
#   ggplot(aes(x = PP.H4)) +
#   geom_histogram(bins = 30, position = "dodge") + 
#   ylim(c(0,550)) + 
#   ggtitle("Distribution of posterior probability of H4 for GEMMAmod's 2,147 gene-peak pairs")
# 
# gemmamod_results |>
#   ggplot(aes(x = PP.H3)) +
#   geom_histogram(bins = 30, position = "dodge") + 
#   ggtitle("Distribution of posterior probability of H3 for GEMMAmod's 2,147 gene-peak pairs")
# 
# gemma_results |>
#   ggplot(aes(x = PP.H4)) +
#   geom_histogram(bins = 30, position = "dodge") + 
#   ylim(c(0,550)) + 
#   ggtitle("Distribution of posterior probability of H4 for GEMMA's 2,120 gene-peak pairs")
# 
# gemma_results |>
#   ggplot(aes(x = PP.H3)) +
#   geom_histogram(bins = 30, position = "dodge") + 
#   ggtitle("Distribution of posterior probability of H3 for GEMMA's 2,120 gene-peak pairs")
# 
# # I want to filter out the gemmamod results that are unique to it
# # and compare those against the GTEx eQTL results
# shared_genepeaks <- comparison[comparison$category == "shared",1:2]
# 
# gemmamod_results_unique <- gemmamod_results_H4H3 %>% 
#   anti_join(shared_genepeaks, by = c("gene_id", "peak_id"))
# 
# gemma_results_unique <- gemma_results_H4H3 %>% 
#   anti_join(shared_genepeaks, by = c("gene_id", "peak_id"))
# 
# # I need to track down the SNPs within these associations
# gemmamod_results_colocSNPs <- gemmamod_eQTL_results_loc %>% 
#   filter(Element_ID %in% gemmamod_results_unique$gene_id) %>% 
#   filter(q_value < 0.05)
# 
# gemma_results_colocSNPs <- gemma_eQTL_results_loc %>% 
#   filter(Element_ID %in% gemma_results_unique$gene_id) %>% 
#   filter(q_value < 0.05)
# 
# # GTEx Pancreas Results
# library(arrow)
# 
# GTEx_pc <- read_parquet("./GTEx_v11/Pancreas.v11.eQTLs.signif_pairs.parquet")
# GTEx_pc_eGenes <- read_tsv("./GTEx_v11/Pancreas.v11.eGenes.txt")
# 
# GTEx_pc_gemmamod_eGenes <- GTEx_pc_eGenes %>% 
#   filter(gene_id %in% gemmamod_results_unique$gene_id) %>% 
#   dplyr::rename(SNP = variant_id) %>% 
#   mutate(SNP = str_remove(SNP, "_b38$"))
# 
# GTEx_pc_gemma_eGenes <- GTEx_pc_eGenes %>% 
#   filter(gene_id %in% gemma_results_unique$gene_id) %>% 
#   dplyr::rename(SNP = variant_id) %>% 
#   mutate(SNP = str_remove(SNP, "_b38$"))
# 
# GTEx_pc_gemmamod_SNPs <- GTEx_pc_gemmamod_geneIDs %>% 
#   filter(SNP %in% gemmamod_results_colocSNPs$SNP)
# 
# GTEx_pc_gemma_SNPs <- GTEx_pc_gemma_geneIDs %>% 
#   filter(SNP %in% gemma_results_colocSNPs$SNP)
# 
# # Perform colocalization between GEMMA/mod and GTEx?
# 
# GTEx_pc_formatted <- GTEx_pc %>% 
#   #filter(phenotype_id %in% gemmamod_eQTL_results_loc$Element_ID) %>% 
#   dplyr::rename(SNP = variant_id,
#                 Element_ID = phenotype_id) %>% 
#   mutate(SNP = str_remove(SNP, "_b38$"))
# GTEx_pc_formatted <- as.data.table(GTEx_pc_formatted)
# 
# # For your colocalized gene-peak pairs, get the contributing SNPs
# # Then check overlap with GTEx significant SNPs
# 
# # Get colocalized genes for each method
# gemmamod_genes <- gemmamod_results_H4H3[, unique(gene_id)]
# gemma_genes    <- gemma_results_H4H3[, unique(gene_id)]
# 
# # Unique to each
# gemmamod_only_genes <- setdiff(gemmamod_genes, gemma_genes)
# gemma_only_genes    <- setdiff(gemma_genes, gemmamod_genes)
# 
# # SNP-level GTEx overlap rate for each set
# calc_gtex_overlap <- function(genes, sig_data, gtex_data) {
#   sig_sub  <- sig_data[Element_ID %in% genes]
#   overlap  <- merge(sig_sub, gtex_data, by = c("Element_ID", "SNP"))
#   n_genes_with_overlap <- uniqueN(overlap$Element_ID)
#   data.table(
#     n_genes        = length(genes),
#     n_with_overlap = n_genes_with_overlap,
#     proportion     = n_genes_with_overlap / length(genes)
#   )
# }
# 
# calc_gtex_overlap(gemmamod_genes, gemmamod_sig, GTEx_pc_formatted)
# calc_gtex_overlap(gemma_genes,    gemma_sig,    GTEx_pc_formatted)
# 
# # Test if the difference is significant
# res_overlap_gemmamod <- calc_gtex_overlap(gemmamod_genes, gemmamod_sig, GTEx_pc_formatted)
# res_overlap_gemma    <- calc_gtex_overlap(gemma_genes,    gemma_sig,    GTEx_pc_formatted)
# 
# prop.test(
#   x = c(res_overlap_gemmamod$n_with_overlap, res_overlap_gemma$n_with_overlap),
#   n = c(res_overlap_gemmamod$n_genes,        res_overlap_gemma$n_genes)
# )
# 
# # Prepare results for QBiC-SELEX analysis
# # Pull all SNPs for the colocalized genes from your significant GEMMAmod results
# gemmamod_coloc_snps <- gemmamod_sig[Element_ID %in% gemmamod_results_H4H3$gene_id]
# gemma_coloc_snps <- gemma_sig[Element_ID %in% gemma_results_H4H3$gene_id]
# 
# # Your SNP IDs look like: chr1_9554508_C_T
# # QBiC needs: chrom, pos, ref, alt
# gemmamod_snp_coords <- gemmamod_coloc_snps[, .(SNP)] |> unique()
# gemma_snp_coords <- gemma_coloc_snps[, .(SNP)] |> unique()
# 
# gemmamod_snp_coords[, c("chrom", "pos", "ref", "alt") := tstrsplit(SNP, "_", fixed = TRUE)]
# gemmamod_snp_coords[, pos := as.integer(pos)]
# 
# gemma_snp_coords[, c("chrom", "pos", "ref", "alt") := tstrsplit(SNP, "_", fixed = TRUE)]
# gemma_snp_coords[, pos := as.integer(pos)]
# 
# # Keep only columns QBiC expects, drop any malformed rows
# gemmamod_qbic_input <- gemmamod_snp_coords[
#   !is.na(pos) & nchar(ref) == 1 & nchar(alt) == 1,  # SNPs only, no indels
#   .(chrom, pos, ref, alt)
# ]
# 
# fwrite(gemmamod_qbic_input, "gemmamod_coloc_variants_for_qbic.csv", col.names = TRUE)
# cat("Variants written:", nrow(gemmamod_qbic_input), "\n")
# 
# gemma_qbic_input <- gemma_snp_coords[
#   !is.na(pos) & nchar(ref) == 1 & nchar(alt) == 1,  # SNPs only, no indels
#   .(chrom, pos, ref, alt)
# ]
# 
# fwrite(gemma_qbic_input, "gemma_coloc_variants_for_qbic.csv", col.names = TRUE)
# cat("Variants written:", nrow(gemma_qbic_input), "\n")
# 
# # ── Read all per-model CSVs from output directory ─────────────────────────────
# gemmamod_qbic_files <- list.files("./QBiC-SELEX/gemmamod_qbic_results/", pattern = "*.csv", full.names = TRUE)
# gemma_qbic_files <- list.files("./QBiC-SELEX/gemma_qbic_results/", pattern = "*.csv", full.names = TRUE)
# 
# library(data.table)
# 
# gemmamod_qbic_files  <- list.files("./QBiC-SELEX/gemmamod_qbic_results/", 
#                           pattern   = "*.csv", 
#                           full.names = TRUE, 
#                           recursive  = TRUE)
# 
# gemmamod_qbic_all <- rbindlist(lapply(gemmamod_qbic_files, fread), fill = TRUE)
# 
# fwrite(gemmamod_qbic_all, "gemmamod_qbic_results_combined.csv")
# gemmamod_qbic_all <- fread(file = "gemmamod_qbic_results_combined.csv")
# 
# gemma_qbic_files  <- list.files("./QBiC-SELEX/gemma_qbic_results/", 
#                                    pattern   = "*.csv", 
#                                    full.names = TRUE, 
#                                    recursive  = TRUE)
# 
# gemma_qbic_all <- rbindlist(lapply(gemma_qbic_files, fread), fill = TRUE)
# 
# fwrite(gemma_qbic_all, "gemma_qbic_results_combined.csv")
# gemma_qbic_all <- fread(file = "gemma_qbic_results_combined.csv")
# 
# # Reconstruct SNP ID in QBiC output to match your format (chr_pos_ref_alt)
# # Adjust column names to match what you actually have
# gemmamod_qbic_all[, SNP := paste(chrom, pos, ref, alt, sep = "_")]
# gemma_qbic_all[, SNP := paste(chrom, pos, ref, alt, sep = "_")]
# 
# # Look at the distribution
# # Combine with method label
# gemmamod_qbic_all[, method := "GEMMAmod"]
# gemma_qbic_all[,    method := "GEMMA"]
# 
# combined_qbic <- rbind(gemmamod_qbic_all, gemma_qbic_all)
# 
# # The number of SNPs with > |5| TF binding change scores between methods
# combined_qbic %>% 
#   group_by(method) %>% 
#   summarize(N_ext = sum(abs(predicted_effect) >= 5))
# 
# nrow(combined_qbic[(combined_qbic$method == "GEMMA" & abs(combined_qbic$predicted_effect) >= 5),])
# 
# prop.test(
#   x = c(nrow(combined_qbic[(combined_qbic$method == "GEMMAmod" & abs(combined_qbic$predicted_effect) >= 5),]), nrow(combined_qbic[(combined_qbic$method == "GEMMA" & abs(combined_qbic$predicted_effect) >= 5),])),
#   n = c(nrow(combined_qbic[(combined_qbic$method == "GEMMAmod"),]),nrow(combined_qbic[(combined_qbic$method == "GEMMA"),]))
# )
# 
# # Overlapping histograms on the same scale
# # T-test first
# ttest_result <- t.test(
#   predicted_effect ~ method,
#   data = combined_qbic
# )
# 
# 
# 
# # Format p-value for plot label
# pval_label <- ifelse(
#   ttest_result$p.value < 0.001,
#   sprintf("t-test p < 0.001 (t = %.2f, df = %.1f)", 
#           ttest_result$statistic, ttest_result$parameter),
#   sprintf("t-test p = %.3f (t = %.2f, df = %.1f)", 
#           ttest_result$p.value, ttest_result$statistic, ttest_result$parameter)
# )
# 
# # Calculate outlier bounds per method
# outliers <- combined_qbic[, {
#   q1  <- quantile(predicted_effect, 0.25, na.rm = TRUE)
#   q3  <- quantile(predicted_effect, 0.75, na.rm = TRUE)
#   iqr <- q3 - q1
#   .SD[predicted_effect < (q1 - 1.5 * iqr) | predicted_effect > (q3 + 1.5 * iqr)]
# }, by = method]
# 
# ggplot(combined_qbic, aes(x = method, y = predicted_effect, color = method)) +
#   geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.4) +
#   geom_jitter(data = outliers, width = 0.2, size = 0.5, alpha = 0.3) +
#   scale_color_manual(values = c("GEMMAmod" = "steelblue", "GEMMA" = "tomato")) +
#   annotate("text", x = 1.5, y = max(combined_qbic$predicted_effect, na.rm = TRUE),
#            label = pval_label, size = 5, hjust = 0.5) +
#   labs(
#     title = "Distribution of QBiC-SELEX Predicted TF Binding Change between GEMMA versions",
#     x     = "Method",
#     y     = "TF binding change: log(mu/wt)",
#     color = "Method"
#   ) +
#   theme_classic() +
#   theme(
#     plot.title   = element_text(hjust = 0.5, size = 14),
#     axis.title   = element_text(size = 13),
#     axis.text    = element_text(size = 12),
#     legend.title = element_text(size = 12),
#     legend.text  = element_text(size = 11)
#   )
# # RUN BELOW SECTION ONCE YOU HAVE STATS
# # ── Filter to significant TF binding changes ──────────────────────────────────
# # p-value threshold is your call; z-score magnitude captures effect size
# qbic_sig <- qbic_results[pvalue < 0.05 & abs(zscore) > 2]
# 
# # ── Reconstruct SNP ID to join back to coloc results ─────────────────────────
# qbic_sig[, SNP := paste(chrom, pos, ref, alt, sep = "_")]
# 
# # ── Join to coloc hits ────────────────────────────────────────────────────────
# coloc_snps[, SNP := paste(chrom, pos, ref, alt, sep = "_")]  # if not already
# 
# qbic_annotated <- merge(
#   qbic_sig,
#   coloc_hits[, .(gene_id, peak_id, PP.H4, PP.H3)],
#   by.x = "SNP",
#   by.y = "SNP",  # adjust if your coloc results use a different SNP key
#   allow.cartesian = TRUE
# )
# 
# fwrite(qbic_annotated, "gemmamod_coloc_qbic_annotated.csv")
# 
# # annotate_coloc_variants.R
# library(data.table)
# library(GenomicRanges)  # BiocManager::install("GenomicRanges")
# library(plyranges)      # BiocManager::install("plyranges") -- tidy GRanges ops
# 
# # ── Load your colocalized variants ────────────────────────────────────────────
# # For GEMMAmod
# gemmamod_snps <- gemmamod_sig[
#   Element_ID %in% gemmamod_results_H4H3$gene_id
# ][, .(SNP)] |> unique()
# gemmamod_snps[, c("chrom", "pos", "ref", "alt") := tstrsplit(SNP, "_")]
# gemmamod_snps[, pos := as.integer(pos)]
# 
# # For GEMMA
# gemma_snps <- gemma_sig[
#   Element_ID %in% gemma_results_H4H3$gene_id
# ][, .(SNP)] |> unique()
# gemma_snps[, c("chrom", "pos", "ref", "alt") := tstrsplit(SNP, "_")]
# gemma_snps[, pos := as.integer(pos)]
# 
# # ── Convert to GRanges (SNPs are single-base, start == end) ───────────────────
# to_granges <- function(dt, method_name) {
#   GRanges(
#     seqnames = dt$chrom,
#     ranges   = IRanges(dt$pos, dt$pos),
#     SNP      = dt$SNP,
#     method   = method_name
#   )
# }
# 
# gr_gemmamod <- to_granges(gemmamod_snps, "GEMMAmod")
# gr_gemma    <- to_granges(gemma_snps,    "GEMMA")
# 
# Panc1Peak <- fread("annotations/wgEncodeRegDnaseUwPanc1Peak.txt.gz")
# colnames(Panc1Peak) <- c("bin", "chrom", "chromStart", "chromEnd", "name", "score",
#                          "strand", "signalValue", "pValue", "qValue", "peak")
# chromhmm_pancreas <- fread(
#   "annotations/roadmap_E087_chromhmm.bed.gz",
#   skip = 1,  # skip the track header line
#   col.names = c("chrom", "start", "end", "state", 
#                 "score", "strand", "thick_start", "thick_end", "color")
# )
# 
# head(chromhmm_pancreas)
# 
# # ── ChromHMM: keep only active regulatory states ──────────────────────────────
# # Quiescent (15), Heterochromatin (9), Repressed Polycomb (13/14) are not functional
# active_states <- c(
#   "1_TssA",      # Active TSS
#   "2_TssAFlnk",  # Flanking Active TSS
#   "3_TxFlnk",    # Transcr. at gene 5' and 3'
#   "6_EnhG",      # Genic Enhancers
#   "7_Enh",       # Enhancers
#   "10_TssBiv",   # Bivalent/Poised TSS -- debatable, remove if you want stricter
#   "12_EnhBiv"    # Bivalent Enhancer   -- debatable, remove if you want stricter
# )
# 
# chromhmm_active <- chromhmm_pancreas[state %in% active_states]
# 
# chromhmm_gr <- GRanges(
#   seqnames = chromhmm_active$chrom,
#   ranges   = IRanges(chromhmm_active$start + 1, chromhmm_active$end),
#   state    = chromhmm_active$state
# )
# 
# # ── ENCODE PANC-1: keep high-confidence peaks only ────────────────────────────
# # signalValue >= 10 removes weak/background peaks
# # pValue >= 10 corresponds to -log10(p) >= 10, i.e. p <= 1e-10
# Panc1Peak_filtered <- Panc1Peak[signalValue >= 10 & pValue >= 10]
# 
# Panc1Peak_gr <- GRanges(
#   seqnames    = Panc1Peak_filtered$chrom,
#   ranges      = IRanges(Panc1Peak_filtered$chromStart + 1, Panc1Peak_filtered$chromEnd),
#   signalValue = Panc1Peak_filtered$signalValue,
#   pValue      = Panc1Peak_filtered$pValue
# )
# 
# # Also break out by category for finer resolution
# chromhmm_promoter <- chromhmm_gr[chromhmm_gr$state %in% c("1_TssA", "2_TssAFlnk")]
# chromhmm_enhancer <- chromhmm_gr[chromhmm_gr$state %in% c("6_EnhG", "7_Enh")]
# chromhmm_bivalent <- chromhmm_gr[chromhmm_gr$state %in% c("10_TssBiv", "12_EnhBiv")]
# 
# annotations <- list(
#   dnase_panc1          = Panc1Peak_gr,
#   chromhmm_promoter    = chromhmm_promoter,
#   chromhmm_enhancer    = chromhmm_enhancer,
#   chromhmm_bivalent    = chromhmm_bivalent
# )
# 
# 
# tally_overlaps <- function(gr_variants, ann_gr, ann_name, method_name) {
#   data.table(
#     annotation    = ann_name,
#     method        = method_name,
#     n_variants    = length(gr_variants),
#     n_overlapping = sum(overlapsAny(gr_variants, ann_gr))
#   )
# }
# 
# tally_results <- rbindlist(c(
#   mapply(tally_overlaps,
#          ann_gr   = annotations,
#          ann_name = names(annotations),
#          MoreArgs = list(gr_variants = gr_gemmamod, method_name = "GEMMAmod"),
#          SIMPLIFY = FALSE
#   ),
#   mapply(tally_overlaps,
#          ann_gr   = annotations,
#          ann_name = names(annotations),
#          MoreArgs = list(gr_variants = gr_gemma, method_name = "GEMMA"),
#          SIMPLIFY = FALSE
#   )
# ))
# 
# tally_results[order(annotation)]
# 
# # Re-run enrichment
# enrichment_results <- rbindlist(mapply(
#   test_enrichment,
#   ann_gr   = annotations,
#   ann_name = names(annotations),
#   MoreArgs = list(gr_gemmamod = gr_gemmamod, gr_gemma = gr_gemma),
#   SIMPLIFY = FALSE
# ))
# 
# enrichment_results[, padj := p.adjust(pvalue, method = "BH")]
# enrichment_results[order(pvalue)]
# 
# 
