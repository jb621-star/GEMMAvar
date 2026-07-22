.libPaths(c("/home/jb621/R/x86_64-pc-linux-gnu-library/4.4","/usr/local/lib/R/site-library","/usr/local/lib/R/library",.libPaths()) ) 
library(parallel)
library(coloc)
library(data.table)
library(tidyverse)

setwd("~/QTL_mapping/")

# Run these next two lines ONLY when you have combinded the two GEMMA version results
# load(file = "PPC_eQTL_SNPs_bothGEMMA.Rda")
# load(file = "PPC_caqtl_SNPs_bothGEMMA.Rda")

load(file = "gemmamod_results_allPairs.Rda")
gemmamod_eQTL_results <- results

load(file = "gemmamod_caqtl_results_allPairs.Rda")
gemmamod_caQTL_results <- results

# Where does _loc come from originally?
# head(gemma_eQTL_results_loc)
# head(gemma_caqtl_results_loc)
# head(gemmamod_eQTL_results_loc)
# head(gemmamod_caqtl_results_loc)
# 
# #save(gemma_eQTL_results_loc, gemmamod_eQTL_results_loc, file = "PPC_eQTL_SNPs_bothGEMMA.Rda")
# #save(gemma_caqtl_results_loc, gemmamod_caqtl_results_loc, file = "PPC_caqtl_SNPs_bothGEMMA.Rda")
# 

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

gemmamod_eQTL_results_loc <- gemmamod_eQTL_results %>%
  dplyr::rename(SNP = SNP_ID)

gemma_eQTL_results_loc <- gemma_eQTL_results %>%
  dplyr::rename(SNP = SNP_ID)

gemmamod_caqtl_results_loc <- gemmamod_caqtl_results %>%
  dplyr::rename(SNP = SNP_ID)

gemma_caqtl_results_loc <- gemma_caqtl_results %>%
  dplyr::rename(SNP = SNP_ID)

# Instead of relying on a set parameter here, need to examine the distribution of p-values
# BiocManager::install("qvalue")
library(qvalue)

# Fit the empirical null using the full p-value distribution
qobj_gemma    <- qvalue(gemma_eQTL_results_loc$p_wald)
qobj_gemmamod <- qvalue(gemmamod_eQTL_results_loc$p_wald)

# Inspect the pi0 estimate (proportion of true nulls)
cat("GEMMA pi0    :", qobj_gemma$pi0, "\n")
cat("GEMMAmod pi0 :", qobj_gemmamod$pi0, "\n")

summary(qobj_gemma)     # should show enrichment near 0 if real signal exists

summary(qobj_gemmamod)     # should show enrichment near 0 if real signal exists

# Add q-values back to your tables
gemma_eQTL_results_loc[,    q_value := qobj_gemma$qvalues]
gemmamod_eQTL_results_loc[, q_value := qobj_gemmamod$qvalues]

# Fit the empirical null using the full p-value distribution
qobj_gemma    <- qvalue(gemma_caqtl_results_loc$p_wald)
qobj_gemmamod <- qvalue(gemmamod_caqtl_results_loc$p_wald)

# Inspect the pi0 estimate (proportion of true nulls)
cat("GEMMA pi0    :", qobj_gemma$pi0, "\n")
cat("GEMMAmod pi0 :", qobj_gemmamod$pi0, "\n")

summary(qobj_gemma)     # should show enrichment near 0 if real signal exists

summary(qobj_gemmamod)     # should show enrichment near 0 if real signal exists

# Add q-values back to your tables
gemma_caqtl_results_loc[,    q_value := qobj_gemma$qvalues]
gemmamod_caqtl_results_loc[, q_value := qobj_gemmamod$qvalues]

# eqtl_snps  <- unique(gemmamod_sig[, .(gene_id = Element_ID, SNP)])
# SNPs tested in caQTL (use union of both methods if they differ)
caqtl_snps <- unique(c(
  gemma_caqtl_results_loc$SNP,
  gemmamod_caqtl_results_loc$SNP
))

# Find overlapping SNPs for each eQTL method
gemma_eQTL_caQTL_SNPs <- gemma_eQTL_results_loc %>%
  dplyr::filter(SNP %in% caqtl_snps) %>%
  dplyr::rename(gene_id = Element_ID)

gemmamod_eQTL_caQTL_SNPs <- gemmamod_eQTL_results_loc %>%
  dplyr::filter(SNP %in% caqtl_snps) %>%
  dplyr::rename(gene_id = Element_ID)

# Step 1: Get the SNPs tested for each gene (already filtered to caQTL SNPs)
gene_snps <- gemmamod_eQTL_caQTL_SNPs[, .(gene_id, SNP)]

# Step 2: Get the SNPs tested for each peak
peak_snps <- gemmamod_caqtl_results_loc[, .(peak_id = Element_ID, SNP)]

# Step 3: Join on SNP — this gives every gene-peak pair that shares at least one SNP
gene_peak_snp <- merge(gene_snps, peak_snps, by = "SNP")

# # Step 4: Check what this looks like
# head(gene_peak_snp)
# cat("Unique genes:", uniqueN(gene_peak_snp$gene_id), "\n")
# cat("Unique peaks:", uniqueN(gene_peak_snp$peak_id), "\n")
# cat("Unique pairs:", nrow(unique(gene_peak_snp[, .(gene_id, peak_id)])), "\n")
# cat("Unique SNPs:", uniqueN(gene_peak_snp$SNP), "\n")

# Collapse to unique pairs
pairs <- unique(gene_peak_snp[, .(gene_id, peak_id)])

library(coloc)
library(data.table)
library(pbapply)

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

save(results, file = "gemmamod_coloc_results_all_pairs.Rda")

# ── Combine GEMMAmod results ───────────────────────────────────────────────────
gemmamod_files <- list.files("coloc_results/",
                             pattern    = "^gemmamod_coloc_results_.*\\.csv$",
                             full.names = TRUE)
cat("GEMMAmod files found:", length(gemmamod_files), "\n")
gemmamod_results <- rbindlist(lapply(gemmamod_files, fread), fill = TRUE)
cat("GEMMAmod results:", nrow(gemmamod_results), "pairs\n")
fwrite(gemmamod_results, "gemmamod_coloc_results_combined.csv")
gemmamod_results <- fread(file = "gemmamod_coloc_results_combined.csv")

# ── Combine GEMMA results ──────────────────────────────────────────────────────
gemma_files <- list.files("coloc_results/",
                          pattern    = "^gemma_coloc_results_.*\\.csv$",
                          full.names = TRUE)
cat("GEMMA files found:", length(gemma_files), "\n")
gemma_results <- rbindlist(lapply(gemma_files, fread), fill = TRUE)
cat("GEMMA results:", nrow(gemma_results), "pairs\n")
fwrite(gemma_results, "gemma_coloc_results_combined.csv")
gemma_results <- fread(file = "gemma_coloc_results_combined.csv")

# Read the results and compare
# load(file = "gemma_coloc_results_all_pairs_051626.Rda")
# gemma_results <- results
# 
# load(file = "gemmamod_coloc_results_all_pairs.Rda")
# gemmamod_results <- results

# # Classify each pair in each result set
classify_coloc <- function(dt, pp4_thresh = 0.75, pp4_3_thresh = 3) {
  dt[, coloc := PP.H4 >= pp4_thresh & (PP.H4 / PP.H3) >= pp4_3_thresh]
  dt
}

gemma_coloc    <- classify_coloc(gemma_results)
gemmamod_coloc <- classify_coloc(gemmamod_results)

# # Merge on gene/peak pair
comparison <- merge(
  gemma_coloc[,    .(gene_id, peak_id, PP.H4.gemma    = PP.H4, PP.H3.gemma    = PP.H3, coloc.gemma    = coloc)],
  gemmamod_coloc[, .(gene_id, peak_id, PP.H4.gemmamod = PP.H4, PP.H3.gemmamod = PP.H3, coloc.gemmamod = coloc)],
  by = c("gene_id", "peak_id")
)

# Four categories
comparison[, category := fcase(
  coloc.gemma &  coloc.gemmamod, "shared",
  coloc.gemma & !coloc.gemmamod, "GEMMA only",
  !coloc.gemma &  coloc.gemmamod, "GEMMAmod only",
  !coloc.gemma & !coloc.gemmamod, "neither"
)]

table(comparison$category)

comparison_shared <- comparison %>% 
  filter(category == "shared") %>% 
  mutate(PP.H4.diff = PP.H4.gemmamod - PP.H4.gemma)

hist(comparison_shared$PP.H4.diff)
abline(v = mean(comparison_shared$PP.H4.diff), col = 'red')

plot(comparison$PP.H4.gemma, comparison$PP.H4.gemmamod,
     main = "Posterior Probability of Gene-Peak Colocalization between GEMMA and GEMMAmod",
     xlab = "PP.H4 for GEMMA",
     ylab = "PP.H4 for GEMMAmod")
abline(a = 0, b = 1, v = 0.75, h = 0.75, col = c('red', 'blue', 'blue'), lty = c(1,2,2), lwd = c(1,2,2))
# What is the distribution of change in confidence?
hist(comparison$PP.H4.gemma - comparison$PP.H4.gemmamod)

sum(comparison_shared$PP.H4.gemmamod >= comparison_shared$PP.H4.gemma)
sum(comparison_shared$PP.H4.gemmamod < comparison_shared$PP.H4.gemma)

# 
gemma_results_H4H3 <- gemma_results %>%
  filter(PP.H4 >= 0.75 & PP.H4 / PP.H3 >= 3)
# 863 / 2120 gene-peak pairs meet the colocalization criteria
length(unique(gemma_results_H4H3$peak_id))
length(unique(gemma_results_H4H3$gene_id))

gemmamod_results_H4H3 <- gemmamod_results %>%
  filter(PP.H4 >= 0.75 & PP.H4 / PP.H3 >= 3)
# 990 / 2147 gene-peak pairs meet the colocalization criteria
length(unique(gemmamod_results_H4H3$peak_id))
length(unique(gemmamod_results_H4H3$gene_id))

gemmamod_results |>
  ggplot(aes(x = PP.H4)) +
  geom_histogram(bins = 30, position = "dodge") +
  #ylim(c(0,550)) +
  ggtitle("Distribution of posterior probability of H4 for GEMMAmod's 2,147 gene-peak pairs")

gemmamod_results |>
  ggplot(aes(x = PP.H3)) +
  geom_histogram(bins = 30, position = "dodge") +
  ggtitle("Distribution of posterior probability of H3 for GEMMAmod's 2,147 gene-peak pairs")

gemma_results |>
  ggplot(aes(x = PP.H4)) +
  geom_histogram(bins = 30, position = "dodge") +
  #ylim(c(0,550)) +
  ggtitle("Distribution of posterior probability of H4 for GEMMA's 2,120 gene-peak pairs")

gemma_results |>
  ggplot(aes(x = PP.H3)) +
  geom_histogram(bins = 30, position = "dodge") +
  ggtitle("Distribution of posterior probability of H3 for GEMMA's 2,120 gene-peak pairs")
# 

head(gemma_caqtl_results_loc)
head(gemmamod_caqtl_results_loc)

# Pull the SNPs that colocalized 
gemmamod_coloc_snps <- gemmamod_caqtl_results_loc[Element_ID %in% gemmamod_results_H4H3$peak_id]
gemma_coloc_snps <- gemma_caqtl_results_loc[Element_ID %in% gemma_results_H4H3$peak_id]

# Your SNP IDs look like: chr1_9554508_C_T
# QBiC needs: chrom, pos, ref, alt
gemmamod_snp_coords <- gemmamod_coloc_snps[, .(SNP)] |> unique()
gemma_snp_coords <- gemma_coloc_snps[, .(SNP)] |> unique()

gemmamod_snp_coords[, c("chrom", "pos", "ref", "alt") := tstrsplit(SNP, "_", fixed = TRUE)]
gemmamod_snp_coords[, pos := as.integer(pos)]

gemma_snp_coords[, c("chrom", "pos", "ref", "alt") := tstrsplit(SNP, "_", fixed = TRUE)]
gemma_snp_coords[, pos := as.integer(pos)]

# Keep only columns QBiC expects, drop any malformed rows
gemmamod_qbic_input <- gemmamod_snp_coords[
  !is.na(pos) & nchar(ref) == 1 & nchar(alt) == 1,  # SNPs only, no indels
  .(chrom, pos, ref, alt)
]

fwrite(gemmamod_qbic_input, "gemmamod_coloc_variants_for_qbic.csv", col.names = TRUE)
cat("Variants written:", nrow(gemmamod_qbic_input), "\n")

gemma_qbic_input <- gemma_snp_coords[
  !is.na(pos) & nchar(ref) == 1 & nchar(alt) == 1,  # SNPs only, no indels
  .(chrom, pos, ref, alt)
]

fwrite(gemma_qbic_input, "gemma_coloc_variants_for_qbic.csv", col.names = TRUE)
cat("Variants written:", nrow(gemma_qbic_input), "\n")

shared_snps <- base::intersect(gemma_snp_coords$SNP, gemmamod_snp_coords$SNP)

# JASPAR TF's intersected with variants

gemmamod_JASPAR_tfs <- fread("./JASPAR/JASPAR2026_gemmamod_variants_overlap.bed") 
gemma_JASPAR_tfs <- fread("./JASPAR/JASPAR2026_gemma_variants_overlap.bed") 

head(gemmamod_coloc_snps)
head(gemmamod_JASPAR_tfs)

# Rename JASPAR columns for clarity
setnames(gemmamod_JASPAR_tfs, 
         c("V1","V2","V3","V4","V5","V6","V7"),
         c("chrom","motif_start","motif_end","motif_ID","score","strand","TF_name"))

# Ensure consistent chromosome column naming for the join
setnames(gemmamod_coloc_snps, "Chromosome", "chrom")

# Non-equi join: SNP position falls within motif window
gemmamod_snp_tf_pairs <- gemmamod_JASPAR_tfs[
  gemmamod_coloc_snps,
  on = .(chrom = chrom, 
         motif_start <= Start, 
         motif_end >= Start),
  nomatch = NULL,
  allow.cartesian = TRUE
][, .(SNP, TF_name, motif_ID, chrom, 
      motif_start, motif_end,    # these now hold the SNP Start position
      Element_ID)]

# Rename JASPAR columns for clarity
setnames(gemma_JASPAR_tfs, 
         c("V1","V2","V3","V4","V5","V6","V7"),
         c("chrom","motif_start","motif_end","motif_ID","score","strand","TF_name"))

# Ensure consistent chromosome column naming for the join
setnames(gemma_coloc_snps, "Chromosome", "chrom")

# Non-equi join: SNP position falls within motif window
gemma_snp_tf_pairs <- gemma_JASPAR_tfs[
  gemma_coloc_snps,
  on = .(chrom = chrom, 
         motif_start <= Start, 
         motif_end >= Start),
  nomatch = NULL,
  allow.cartesian = TRUE
][, .(SNP, TF_name, motif_ID, chrom, 
      motif_start, motif_end,    # these now hold the SNP Start position
      Element_ID)]



# Roadmap for the rest of this analysis is to export each set of SNPs against
gemmamod_tfs <- unique(gemmamod_snp_tf_pairs$TF_name)
gemma_tfs <- unique(gemma_snp_tf_pairs$TF_name)

# Get unique TF names from your pairs table
gemmamod_tf_names <- unique(gemmamod_snp_tf_pairs$TF_name)
gemma_tf_names <- unique(gemma_snp_tf_pairs$TF_name)

# Filter out any TF names containing "::" (heterodimer conventions)
gemmamod_tf_names_clean <- gemmamod_tf_names[!grepl("::", gemmamod_tf_names, fixed = TRUE)]
gemma_tf_names_clean <- gemma_tf_names[!grepl("::", gemma_tf_names, fixed = TRUE)]

for (tf in gemmamod_tf_names_clean) {
  #tf <- "NHLH2"
  tf_upper <- toupper(tf)
  
  gemmamod_tf_snps <- unique(gemmamod_snp_tf_pairs %>% 
    filter(TF_name == tf) %>% 
    select(SNP))
  
  gemmamod_qbic_input <- gemmamod_tf_snps %>% 
    separate(SNP, into = c("chrom", "pos", "ref", "alt"))
  
  fwrite(gemmamod_qbic_input, paste0("./QBiC-SELEX/gemmamod_coloc_variants_for_", tf_upper,".csv"), col.names = TRUE)
}

for (tf in gemma_tf_names_clean) {
  tf_upper <- toupper(tf)
  
  gemma_tf_snps <- unique(gemma_snp_tf_pairs %>% 
                               filter(TF_name == tf) %>% 
                               select(SNP))
  
  gemma_qbic_input <- gemma_tf_snps %>% 
    separate(SNP, into = c("chrom", "pos", "ref", "alt"))
  
  fwrite(gemma_qbic_input, paste0("./QBiC-SELEX/gemma_coloc_variants_for_", tf_upper,".csv"), col.names = TRUE)
  
  
}

# # ── Read all per-model CSVs from output directory ─────────────────────────────
gemmamod_qbic_files  <- list.files("./QBiC-SELEX/gemmamod_qbic_results_v2/",
                          pattern   = "*.csv",
                          full.names = TRUE,
                          recursive  = TRUE)

gemmamod_qbic_all <- rbindlist(lapply(gemmamod_qbic_files, fread), fill = TRUE)

gemmamod_qbic_all <- gemmamod_qbic_all %>% 
  mutate(SNP = paste0(chrom, "_", pos, "_", ref,"_", alt))

fwrite(gemmamod_qbic_all, "gemmamod_qbic_results_v2_combined.csv")
gemmamod_qbic_all <- fread(file = "gemmamod_qbic_results_v2_combined.csv")

gemma_qbic_files  <- list.files("./QBiC-SELEX/gemma_qbic_results_v2/",
                                   pattern   = "*.csv",
                                   full.names = TRUE,
                                   recursive  = TRUE)

gemma_qbic_all <- rbindlist(lapply(gemma_qbic_files, fread), fill = TRUE)

gemma_qbic_all <- gemma_qbic_all %>% 
  mutate(SNP = paste0(chrom, "_", pos, "_", ref,"_", alt))

fwrite(gemma_qbic_all, "gemma_qbic_results_v2_combined.csv")
gemma_qbic_all <- fread(file = "gemma_qbic_results_v2_combined.csv")

# RUN BELOW SECTION ONCE YOU HAVE STATS
# gemma_qbic_all has 3,200,967 observations
# gemmamod_qbic_all has 3,946,734 observations
# This nummber for GEMMAmod is consistent with 729 more SNPs X 1023 TF motifs
gemma_qbic_all[, SNP := paste(chrom, pos, ref, alt, sep = "_")]
gemma_qbic_all_u <- gemma_qbic_all %>% 
  filter(!(SNP %in% shared_snps))

gemmamod_qbic_all[, SNP := paste(chrom, pos, ref, alt, sep = "_")]
gemmamod_qbic_all_u <- gemmamod_qbic_all %>% 
  filter(!(SNP %in% shared_snps))

# ── Filter to significant TF binding changes ──────────────────────────────────
# p-value threshold is your call; z-score magnitude captures effect size
gemma_qbic_sig <- gemma_qbic_all[p_value < 0.05 & abs(z_score) > 2]
gemmamod_qbic_sig <- gemmamod_qbic_all[p_value < 0.05 & abs(z_score) > 2]

# ── Reconstruct SNP ID to join back to coloc results ─────────────────────────
gemma_qbic_sig[, SNP := paste(chrom, pos, ref, alt, sep = "_")]
gemma_qbic_sig_u <- gemma_qbic_sig %>% 
  filter(!(SNP %in% shared_snps))
hist(gemma_qbic_sig_u$predicted_effect)
sum(abs(gemma_qbic_sig_u$predicted_effect) >= 5)  
# 283,392 / 435,798 significant SNP/TF interactions unique to GEMMA (65.0%)
# 96 / 283,392 (0.03%) had a |predicted effect| >= 5 

gemmamod_qbic_sig[, SNP := paste(chrom, pos, ref, alt, sep = "_")]
gemmamod_qbic_sig_u <- gemmamod_qbic_sig %>% 
  filter(!(SNP %in% shared_snps))
hist(gemmamod_qbic_sig_u$predicted_effect)
sum(abs(gemmamod_qbic_sig_u$predicted_effect) >= 5)  
# 757,615 / 1,181,565 significant SNP/TF interactions unique to GEMMAmod (64.1%)
# 275 / 757,615 (0.04%) had a |predicted effect| >= 5 

# Is gemmamod_qbic_sig = gemma + more?
sum(gemma_qbic_sig$SNP %in% gemmamod_qbic_sig$SNP)
# gemma and gemmamod share 1,777,470 entries

# ── Join to coloc hits ────────────────────────────────────────────────────────
# ── Tag each dataset before merging ──────────────────────────────────────────
gemma_qbic_sig[, source := "GEMMA"]
gemmamod_qbic_sig[, source := "GEMMAmod"]

# ── Create a SNP_TF key for matching across both datasets ────────────────────
gemma_qbic_sig[,    SNP_TF := paste(SNP, model, sep = "_")]
gemmamod_qbic_sig[, SNP_TF := paste(SNP, model, sep = "_")]

# ── Find overlapping SNP/TF pairs ────────────────────────────────────────────
shared_snp_tf <- intersect(gemma_qbic_sig$SNP_TF, gemmamod_qbic_sig$SNP_TF)

# ── Tag shared pairs in both datasets ────────────────────────────────────────
gemma_qbic_sig[,    source := ifelse(SNP_TF %in% shared_snp_tf, "both", "GEMMA")]
gemmamod_qbic_sig[, source := ifelse(SNP_TF %in% shared_snp_tf, "both", "GEMMAmod")]

# ── Merge, keeping only one copy of shared rows ──────────────────────────────
# Use GEMMAmod's values for shared rows (or GEMMA's — your call)
qbic_merged <- rbind(
  gemma_qbic_sig[source == "GEMMA"],        # unique to GEMMA
  gemmamod_qbic_sig[source == "GEMMAmod"],  # unique to GEMMAmod
  gemmamod_qbic_sig[source == "both"]       # shared (using GEMMAmod values)
)

# ── Sanity check ─────────────────────────────────────────────────────────────
table(qbic_merged$source)

# fwrite(qbic_annotated, "gemmamod_coloc_qbic_annotated.csv")
# 

# Do all the JASPR TFs overlap with ATAC peaks?
head(gemmamod_JASPAR_tfs)
gemmamod_JASPAR_tfs %>% 
  filter(chrom == "chr10" & motif_start <= 101013254 & motif_end >= 101013254)

# Colocalizing instead of just correlating
compute_locus_coloc <- function(results_H4H3, eQTL_SNPs, caqtl_SNPs, qbic_all, method_name, N_samples) {
  # results_H4H3 <- gemma_results_H4H3
  # eQTL_SNPs <- gemma_eQTL_caQTL_SNPs
  # caqtl_SNPs <- gemma_caqtl_results_loc
  # qbic_all <- gemma_qbic_all_tf
  # method_name <- "GEMMA"
  # N_samples <- 105
  coloc_results <- list()
  
  for (i in 1:nrow(results_H4H3)) {
    #i <- 1
    gene <- results_H4H3$gene_id[i]
    peak <- results_H4H3$peak_id[i]
    
    # Get caQTL summary stats at this locus
    locus_caqtl <- caqtl_SNPs %>%
      filter(Element_ID == peak) %>%
      group_by(SNP) #%>%
      #slice(which.max(-log10(q_value))) %>%
      #ungroup()
    
    # Get eQTL summary stats at this locus
    locus_eqtl <- eQTL_SNPs %>%
      filter(gene_id == gene) %>%
      filter(SNP %in% locus_caqtl$SNP)
    
    # Get QBiC scores
    locus_snps <- intersect(locus_eqtl$SNP, locus_caqtl$SNP)
    locus_qbic <- qbic_all %>%
      filter(SNP %in% locus_snps)
    
    if (nrow(locus_qbic) > 0) {
    
    # ── eQTL vs QBiC coloc ───────────────────────────────────────────────────
    eqtl_qbic <- inner_join(locus_eqtl, locus_qbic, by = "SNP")
    # QBiC provides: predicted_effect (cTβ̂) and z-score
    # From these you can derive a pseudo-SE:
    # SE = predicted_effect / z_score
    eqtl_qbic <- eqtl_qbic %>% 
      mutate(pseudo_se = predicted_effect / z_score)
      
      d1_eqtl <- list(
        beta    = eqtl_qbic$beta,
        varbeta = eqtl_qbic$se_beta^2,
        N       = N_samples,
        sdY     = sd(t(edgeR_expr_ordered[gene, ])),
        type    = "quant",
        snp     = eqtl_qbic$SNP
      )
      
      d2_qbic <- list(
        beta    = eqtl_qbic$predicted_effect,
        varbeta = eqtl_qbic$pseudo_se^2,
        N       = N_samples,
        sdY     = 1,
        type    = "quant",
        snp     = eqtl_qbic$SNP
      )
      
      eqtl_coloc <- tryCatch({
        coloc.abf(d1_eqtl, d2_qbic)$summary
      }, error = function(e) NULL)
      
    
    # ── caQTL vs QBiC coloc ──────────────────────────────────────────────────
    caqtl_qbic <- inner_join(locus_caqtl, locus_qbic, by = "SNP")
    # QBiC provides: predicted_effect (cTβ̂) and z-score
    # From these you can derive a pseudo-SE:
    # SE = predicted_effect / z_score
    caqtl_qbic <- caqtl_qbic %>% 
      mutate(pseudo_se = predicted_effect / z_score)
    
    #if (nrow(caqtl_qbic) >= 4) {
      
      d1_caqtl <- list(
        beta    = caqtl_qbic$beta,
        varbeta = caqtl_qbic$se_beta^2,
        N       = N_samples,
        sdY     = sd(t(atac_expr_ordered[peak,  ])),
        type    = "quant",
        snp     = caqtl_qbic$SNP
      )

      d2_qbic_ca <- list(
        beta    = caqtl_qbic$predicted_effect,
        varbeta = caqtl_qbic$pseudo_se^2,
        N       = N_samples,
        sdY     = 1,
        type    = "quant",
        snp     = caqtl_qbic$SNP
      )
      
      caqtl_coloc <- tryCatch({
        coloc.abf(d1_caqtl, d2_qbic_ca)$summary
      }, error = function(e) NULL)
      
    } else {
     eqtl_coloc <- NULL
     caqtl_coloc <- NULL
    }
    
    # ── Store results ─────────────────────────────────────────────────────────
    coloc_results[[i]] <- data.frame(
      gene_id        = gene,
      peak_id        = peak,
      method         = method_name,
      eqtl_nsnps     = if (!is.null(eqtl_coloc))  nrow(eqtl_qbic)           else NA,
      eqtl_PP.H0     = if (!is.null(eqtl_coloc))  eqtl_coloc["PP.H0.abf"]   else NA,
      eqtl_PP.H3     = if (!is.null(eqtl_coloc))  eqtl_coloc["PP.H3.abf"]   else NA,
      eqtl_PP.H4     = if (!is.null(eqtl_coloc))  eqtl_coloc["PP.H4.abf"]   else NA,
      caqtl_nsnps    = if (!is.null(caqtl_coloc)) nrow(caqtl_qbic)           else NA,
      caqtl_PP.H0    = if (!is.null(caqtl_coloc)) caqtl_coloc["PP.H0.abf"]   else NA,
      caqtl_PP.H3    = if (!is.null(caqtl_coloc)) caqtl_coloc["PP.H3.abf"]   else NA,
      caqtl_PP.H4    = if (!is.null(caqtl_coloc)) caqtl_coloc["PP.H4.abf"]   else NA
    )
  }
  
  bind_rows(coloc_results)
}

PPC_TFs <- c("PDX1", "NKX6-1", "PTF1A", "SOX9", 
             "NGN3", "NKX2-2", "GLIS3", "HES1",
             "FOXA2", "FOXA3", "ONECUT3", "GCK",
             "NEUROG3", "PAX6", "HNF1A", "HNF4A",
             "HNF4G", "PROX1", "INS", "INSM1",
             "NEUROD1", "NRSA2", "HLXB9", "HNF1B", "GATA6")

# Initialize storage OUTSIDE the loop
gemma_tf_coloc_results    <- list()
gemmamod_tf_coloc_results <- list()

for (tf in PPC_TFs) {
  
  gemma_qbic_all_tf <- gemma_qbic_all %>% 
    filter(grepl(tf, model))
  
  gemmamod_qbic_all_tf <- gemmamod_qbic_all %>% 
    filter(grepl(tf, model))
  
  # Skip if no QBiC results exist for this TF
  if (nrow(gemma_qbic_all_tf) == 0 & nrow(gemmamod_qbic_all_tf) == 0) {
    message("No QBiC results found for TF: ", tf, " — skipping.")
    next
  }
  
  gemma_qbic_coloc_tf <- compute_locus_coloc(
    gemma_results_H4H3,
    gemma_eQTL_caQTL_SNPs,
    gemma_caqtl_results_loc,
    gemma_qbic_all_tf,
    "GEMMA",
    N_samples = 105
  )
  
  gemmamod_qbic_coloc_tf <- compute_locus_coloc(
    gemmamod_results_H4H3,
    gemmamod_eQTL_caQTL_SNPs,
    gemmamod_caqtl_results_loc,
    gemmamod_qbic_all_tf,
    "GEMMAmod",
    N_samples = 105
  )
  
  # Tag results with TF name and store
  gemma_tf_coloc_results[[tf]]    <- gemma_qbic_coloc_tf    %>% mutate(TF = tf)
  gemmamod_tf_coloc_results[[tf]] <- gemmamod_qbic_coloc_tf %>% mutate(TF = tf)
  
  message("Completed TF: ", tf)
}

# Combine into single dataframes
gemma_tf_coloc_all    <- bind_rows(gemma_tf_coloc_results)
gemmamod_tf_coloc_all <- bind_rows(gemmamod_tf_coloc_results)

# Inspect
print(gemma_tf_coloc_all)
print(gemmamod_tf_coloc_all)

gemma_tf_coloc_all_filt <- gemma_tf_coloc_all %>% 
  filter(eqtl_nsnps >= 2)

gemmamod_tf_coloc_all_filt <- gemmamod_tf_coloc_all %>% 
  filter(eqtl_nsnps >= 2)

save(gemma_tf_coloc_all_filt, file = "gemma_tf_coloc_all_0610.Rda")
load(file = "gemma_tf_coloc_all_0610.Rda")

save(gemmamod_tf_coloc_all_filt, file = "gemmamod_tf_coloc_all_0610.Rda")
load(file = "gemmamod_tf_coloc_all_0610.Rda")

# ── Run for both methods ──────────────────────────────────────────────────────
# N_samples = number of individuals in your QTL analysis
# gemma_qbic_coloc <- compute_locus_coloc(
#   gemma_results_H4H3,
#   gemma_eQTL_caQTL_SNPs,
#   gemma_caqtl_results_loc,
#   gemma_qbic_all,
#   "GEMMA",
#   N_samples = 105
# )
gemma_qbic_coloc <- gemma_tf_coloc_all_filt %>% 
  mutate(eqtl_coloc = if_else((eqtl_PP.H4 >= 0.75 & eqtl_PP.H4 / eqtl_PP.H3 >= 3), TRUE, FALSE),
         caqtl_coloc = if_else((caqtl_PP.H4 >= 0.75 & caqtl_PP.H4 / caqtl_PP.H3 >= 3), TRUE, FALSE))

# gemmamod_qbic_coloc <- compute_locus_coloc(
#   gemmamod_results_H4H3,
#   gemmamod_eQTL_caQTL_SNPs,
#   gemmamod_caqtl_results_loc,
#   gemmamod_qbic_all,
#   "GEMMAmod",
#   N_samples = 105
# )
gemmamod_qbic_coloc <- gemmamod_tf_coloc_all_filt %>% 
  mutate(eqtl_coloc = if_else((eqtl_PP.H4 >= 0.75 & eqtl_PP.H4 / eqtl_PP.H3 >= 3), TRUE, FALSE),
         caqtl_coloc = if_else((caqtl_PP.H4 >= 0.75 & caqtl_PP.H4 / caqtl_PP.H3 >= 3), TRUE, FALSE))


qbic_tf_coloc_all <- bind_rows(gemma_tf_coloc_all_filt, gemmamod_tf_coloc_all_filt)

qbic_tf_coloc_all <- qbic_tf_coloc_all %>% 
  mutate(eqtl_coloc = if_else((eqtl_PP.H4 >= 0.75 & eqtl_PP.H4 / eqtl_PP.H3 >= 3), TRUE, FALSE),
         caqtl_coloc = if_else((caqtl_PP.H4 >= 0.75 & caqtl_PP.H4 / caqtl_PP.H3 >= 3), TRUE, FALSE))

qbic_tf_coloc_all_table <- qbic_tf_coloc_all %>% 
  group_by(method, TF) %>% 
  summarise(eqtl_coloc_rate = sum(eqtl_coloc),
            caqtl_coloc_rate = sum(caqtl_coloc))

table(qbic_tf_coloc_all$method, qbic_tf_coloc_all$eqtl_coloc)
table(qbic_tf_coloc_all$method, qbic_tf_coloc_all$caqtl_coloc)


# Filter to the shared gene/peak pairs, and see what the colocalization probabilities are
save(qbic_tf_coloc_all, file = "qbic_tf_coloc_all_v2_0610.Rda")
load(file = "qbic_tf_coloc_all_v2_0610.Rda")

gemma_qbic_coloc_pairs <- paste0(gemma_qbic_coloc$gene_id, "_", gemma_qbic_coloc$peak_id)
gemmamod_qbic_coloc_pairs <- paste0(gemmamod_qbic_coloc$gene_id, "_", gemmamod_qbic_coloc$peak_id)

shared_qbic_coloc_pairs <- intersect(gemma_qbic_coloc_pairs, gemmamod_qbic_coloc_pairs)

qbic_coloc_shared <- qbic_tf_coloc_all %>% 
  mutate(pair_key = paste0(gene_id, "_", peak_id)) %>% 
  filter(pair_key %in% shared_qbic_coloc_pairs) %>% 
  mutate(eqtl_coloc = if_else((eqtl_PP.H4 >= 0.75 & eqtl_PP.H4 / eqtl_PP.H3 >= 3), TRUE, FALSE),
         caqtl_coloc = if_else((caqtl_PP.H4 >= 0.75 & caqtl_PP.H4 / caqtl_PP.H3 >= 3), TRUE, FALSE))


qbic_coloc_shared %>% 
  group_by(method) %>% 
  summarise(eqtl_coloc_rate = sum(eqtl_coloc) / length(shared_qbic_coloc_pairs),
            caqtl_coloc_rate = sum(caqtl_coloc) / length(shared_qbic_coloc_pairs))

# Build a unified contingency table for McNemar's Test
# GEMMA ver
head(gemma_coloc)
head(gemma_qbic_coloc)

gemma_qbic_coloc_TF <- gemma_qbic_coloc %>% 
  filter(TF == "SOX9")

library(dplyr)

# Start with ALL stage-1 tested pairs
gemma_full <- gemma_coloc %>%
  select(gene_id, peak_id) %>%
  distinct() %>%
  
  # Add stage-2 results where available
  left_join(
    gemma_qbic_coloc_TF %>%
      select(
        gene_id,
        peak_id,
        eqtl_coloc,
        caqtl_coloc
      ),
    by = c("gene_id", "peak_id")
  ) %>%
  
  # Replace missing stage-2 results with FALSE
  mutate(
    eqtl_coloc  = ifelse(is.na(eqtl_coloc), FALSE, eqtl_coloc),
    caqtl_coloc = ifelse(is.na(caqtl_coloc), FALSE, caqtl_coloc)
  )

gemma_full %>% 
  summarise(eqtl_coloc_rate = sum(eqtl_coloc),
            caqtl_coloc_rate = sum(caqtl_coloc))

# GEMMAmod ver

gemmamod_qbic_coloc_TF <- gemmamod_qbic_coloc %>% 
  filter(TF == "SOX9")

# Start with ALL stage-1 tested pairs
gemmamod_full <- gemmamod_coloc %>%
  select(gene_id, peak_id) %>%
  distinct() %>%
  
  # Add stage-2 results where available
  left_join(
    gemmamod_qbic_coloc_TF %>%
      select(
        gene_id,
        peak_id,
        eqtl_coloc,
        caqtl_coloc
      ),
    by = c("gene_id", "peak_id")
  ) %>%
  
  # Replace missing stage-2 results with FALSE
  mutate(
    eqtl_coloc  = ifelse(is.na(eqtl_coloc), FALSE, eqtl_coloc),
    caqtl_coloc = ifelse(is.na(caqtl_coloc), FALSE, caqtl_coloc)
  )

gemmamod_full %>% 
  summarise(eqtl_coloc_rate = sum(eqtl_coloc),
            caqtl_coloc_rate = sum(caqtl_coloc))


eqtl_paired <- inner_join(
  gemma_full %>%
    select(gene_id, peak_id, gemma = eqtl_coloc),
  
  gemmamod_full %>%
    select(gene_id, peak_id, gemmamod = eqtl_coloc),
  
  by = c("gene_id", "peak_id")
)

eqtl_tab <- table(
  eqtl_paired$gemma,
  eqtl_paired$gemmamod
)

eqtl_tab

mcnemar.test(eqtl_tab)

caqtl_paired <- inner_join(
  gemma_full %>%
    select(gene_id, peak_id, gemma = caqtl_coloc),
  
  gemmamod_full %>%
    select(gene_id, peak_id, gemmamod = caqtl_coloc),
  
  by = c("gene_id", "peak_id")
)

caqtl_tab <- table(
  caqtl_paired$gemma,
  caqtl_paired$gemmamod
)

caqtl_tab

mcnemar.test(caqtl_tab)


# Genome tracks for the same SNPs and TF factors that iPSCORE outlined
gene_of_interest <- "ENSG00000153814.13"
gene_of_interest <- "ENSG00000079689.15"

gene_of_interest <- "ENSG00000109445.11"
gene_of_interest <- "ENSG00000042980.14"
gene_of_interest <- "ENSG00000013725.14"
gene_of_interest <- "ENSG00000250305.9"

peak_of_interest <- "ppc_atac_peak_244298"
peak_of_interest <- "ppc_atac_peak_227797"

peak_of_interest <- "ppc_atac_peak_203908"
peak_of_interest <- "ppc_atac_peak_258526"
peak_of_interest <- "ppc_atac_peak_47038"
peak_of_interest <- "ppc_atac_peak_257307"

# Need to find a gene/peak pair with GEMMAmod eqtl/caqtl coloc but not GEMMA
# Find gene_id/peak_id pairs where GEMMA has both FALSE
# but GEMMAmod has both TRUE

# Check eqtl only
gemma_eqtl_false <- qbic_tf_coloc_all %>%
  filter(method == "GEMMA", eqtl_coloc == FALSE) %>%
  select(gene_id, peak_id, TF)

gemma_eqtl_true <- qbic_tf_coloc_all %>%
  filter(method == "GEMMA", eqtl_coloc == TRUE) %>%
  select(gene_id, peak_id, TF)

gemmamod_eqtl_true <- qbic_tf_coloc_all %>%
  filter(method == "GEMMAmod", eqtl_coloc == TRUE) %>%
  select(gene_id, peak_id, TF)

gemmamod_eqtl_false <- qbic_tf_coloc_all %>%
  filter(method == "GEMMAmod", eqtl_coloc == FALSE) %>%
  select(gene_id, peak_id, TF)

eqtl_switched_gemmaT <- inner_join(gemma_eqtl_true, gemmamod_eqtl_false, by = c("gene_id", "peak_id", "TF"))
eqtl_switched_gemmaF <- inner_join(gemma_eqtl_false, gemmamod_eqtl_true, by = c("gene_id", "peak_id", "TF"))
cat("eQTL only switched:", nrow(eqtl_switched), "\n")

eqtl_con_T <- inner_join(gemma_eqtl_true, gemmamod_eqtl_true, by = c("gene_id", "peak_id", "TF"))
eqtl_con_F <- inner_join(gemma_eqtl_false, gemmamod_eqtl_false, by = c("gene_id", "peak_id", "TF"))


# Check caqtl only
gemma_caqtl_false <- qbic_tf_coloc_all %>%
  filter(method == "GEMMA", caqtl_coloc == FALSE) %>%
  select(gene_id, peak_id, TF)

gemma_caqtl_true <- qbic_tf_coloc_all %>%
  filter(method == "GEMMA", caqtl_coloc == TRUE) %>%
  select(gene_id, peak_id, TF)

gemmamod_caqtl_true <- qbic_tf_coloc_all %>%
  filter(method == "GEMMAmod", caqtl_coloc == TRUE) %>%
  select(gene_id, peak_id, TF)

gemmamod_caqtl_false <- qbic_tf_coloc_all %>%
  filter(method == "GEMMAmod", caqtl_coloc == FALSE) %>%
  select(gene_id, peak_id, TF)

caqtl_switched_gemmaF <- inner_join(gemma_caqtl_false, gemmamod_caqtl_true, by = c("gene_id", "peak_id", "TF"))
caqtl_switched_gemmaT <- inner_join(gemma_caqtl_true, gemmamod_caqtl_false, by = c("gene_id", "peak_id", "TF"))
cat("caQTL only switched:", nrow(caqtl_switched), "\n")

caqtl_con_T <- inner_join(gemma_caqtl_true, gemmamod_caqtl_true, by = c("gene_id", "peak_id", "TF"))
caqtl_con_F <- inner_join(gemma_caqtl_false, gemmamod_caqtl_false, by = c("gene_id", "peak_id", "TF"))



# Show results
eqtl_switched
caqtl_switched

# GEMMA
gemma_locus_snps_g_eQTL <- gemma_eQTL_caQTL_SNPs %>% 
  filter(gene_id == gene_of_interest) %>% 
  dplyr::arrange(Start)

gemma_locus_snps_g_caQTL <- gemma_caqtl_results_loc %>% 
  filter(Element_ID == peak_of_interest) %>% 
  dplyr::arrange(Start)

gemma_locus_snps_g_shared <- intersect(
  gemma_locus_snps_g_eQTL$SNP,
  gemma_locus_snps_g_caQTL$SNP
)

# Filter eQTL results to only those SNPs for a fair comparison
gemma_locus_snps_g_eQTL <- gemma_locus_snps_g_eQTL %>% 
  filter(SNP %in% gemma_locus_snps_g_shared)

gemma_locus_snps_g_caQTL <- gemma_locus_snps_g_caQTL %>% 
  filter(SNP %in% gemma_locus_snps_g_shared)

# GEMMAmod
gemmamod_locus_snps_g_eQTL <- gemmamod_eQTL_caQTL_SNPs %>% 
  filter(gene_id == gene_of_interest) %>% 
  dplyr::arrange(Start)

gemmamod_locus_snps_g_caQTL <- gemmamod_caqtl_results_loc %>% 
  filter(Element_ID == peak_of_interest) %>% 
  dplyr::arrange(Start)

gemmamod_locus_snps_g_shared <- intersect(
  gemmamod_locus_snps_g_eQTL$SNP,
  gemmamod_locus_snps_g_caQTL$SNP
)

# Filter eQTL results to only those SNPs for a fair comparison
gemmamod_locus_snps_g_eQTL <- gemmamod_locus_snps_g_eQTL %>% 
  filter(SNP %in% gemmamod_locus_snps_g_shared)

gemmamod_locus_snps_g_caQTL <- gemmamod_locus_snps_g_caQTL %>% 
  filter(SNP %in% gemmamod_locus_snps_g_shared)

locus_snps_g_eQTL <- data.frame(SNP = gemma_locus_snps_g_eQTL$SNP,
                                 Chromosome = gemma_locus_snps_g_eQTL$Chromosome,
                                 Start = gemma_locus_snps_g_eQTL$Start,
                                 beta_def = gemma_locus_snps_g_eQTL$beta,
                                 beta_se_def = gemma_locus_snps_g_eQTL$se_beta,
                                 q_value_def = gemma_locus_snps_g_eQTL$q_value,
                                log10q_def = -log10(gemma_locus_snps_g_eQTL$q_value),
                                 beta_mod = gemmamod_locus_snps_g_eQTL$beta,
                                 beta_se_mod = gemmamod_locus_snps_g_eQTL$se_beta,
                                 q_value_mod = gemmamod_locus_snps_g_eQTL$q_value,
                                log10q_mod = -log10(gemmamod_locus_snps_g_eQTL$q_value))

locus_snps_g_caQTL <- data.frame(SNP = gemma_locus_snps_g_caQTL$SNP,
                                 Chromosome = gemma_locus_snps_g_caQTL$Chromosome,
                                 Start = gemma_locus_snps_g_caQTL$Start,
                                 beta_def = gemma_locus_snps_g_caQTL$beta,
                                 beta_se_def = gemma_locus_snps_g_caQTL$se_beta,
                                 q_value_def = gemma_locus_snps_g_caQTL$q_value,
                                 log10q_def = -log10(gemma_locus_snps_g_caQTL$q_value),
                                 beta_mod = gemmamod_locus_snps_g_caQTL$beta,
                                 beta_se_mod = gemmamod_locus_snps_g_caQTL$se_beta,
                                 q_value_mod = gemmamod_locus_snps_g_caQTL$q_value,
                                 log10q_mod = -log10(gemmamod_locus_snps_g_caQTL$q_value))

# --- Reshape eQTL data to long format ---
locus_snps_g_eQTL_long <- locus_snps_g_eQTL %>%
  pivot_longer(
    cols = c(log10q_def, log10q_mod),
    names_to = "method",
    values_to = "log10q"
  ) %>%
  mutate(method = recode(method,
                         "log10q_def" = "GEMMA",
                         "log10q_mod" = "GEMMAmod"))

# --- Reshape caQTL data to long format ---
locus_snps_g_caQTL_long <- locus_snps_g_caQTL %>%
  pivot_longer(
    cols = c(log10q_def, log10q_mod),
    names_to = "method",
    values_to = "log10q"
  ) %>%
  mutate(method = recode(method,
                         "log10q_def" = "GEMMA",
                         "log10q_mod" = "GEMMAmod"))

# --- eQTL plot ---
eQTL_results <- ggplot(data = locus_snps_g_eQTL_long,
                       aes(x = Start, y = log10q, color = method)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("GEMMA" = "red", "GEMMAmod" = "blue"),
                     name = "Method") +
  ggtitle(paste0("GEMMA vs GEMMAmod eQTL results for ",gene_of_interest)) +
  ylab("-log10(FDR adj. p-value)") +
  scale_x_continuous("Genomic Position (Mb)", labels = scales::label_number(scale = 1e-6)) +
  facet_grid(. ~ Chromosome, scales = "free_x") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

eQTL_results

# --- caQTL plot ---
caQTL_results <- ggplot(data = locus_snps_g_caQTL_long,
                        aes(x = Start, y = log10q, color = method)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("GEMMA" = "red", "GEMMAmod" = "blue"),
                     name = "Method") +
  ggtitle(paste0("GEMMA vs GEMMAmod caQTL results for ",peak_of_interest)) +
  ylab("-log10(FDR adj. p-value)") +
  scale_x_continuous("Genomic Position (Mb)", labels = scales::label_number(scale = 1e-6)) +
  facet_grid(. ~ Chromosome, scales = "free_x") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

caQTL_results

gemma_eQTL_caQTL_SNPs %>% 
  filter(SNP == "chr6_25660580_C_T" & gene_id == "ENSG00000079689.15")

gemmamod_eQTL_caQTL_SNPs %>% 
  filter(SNP == "chr6_25660580_C_T" & gene_id == "ENSG00000079689.15")

gemma_eQTL_caQTL_SNPs %>% 
  filter(SNP == "chr6_25661140_A_C" & gene_id == "ENSG00000079689.15")

gemmamod_eQTL_caQTL_SNPs %>% 
  filter(SNP == "chr6_25661140_A_C" & gene_id == "ENSG00000079689.15")


# Pull the predicted variant effects for these variants for PDX1 and NKX6-1
gemmamod_qbic_all %>% 
  filter(SNP %in% gemmamod_locus_snps_g_caQTL$SNP) %>% 
  filter(grepl("NKX6-1", model))

gemmamod_qbic_all %>% 
  filter(SNP %in% gemmamod_locus_snps_g_caQTL$SNP) %>% 
  filter(grepl("PTF1A", model))

# PDX1, NKX6-1 TF model results for these SNPs
gemmamod_qbic_all %>% 
  filter(SNP %in% gemmamod_locus_snps_g_caQTL$SNP) %>% 
  filter(grepl("NKX6-1", model))

gemmamod_qbic_all %>% 
  filter(SNP %in% gemmamod_locus_snps_g_caQTL$SNP) %>% 
  filter(grepl("PDX1", model))

gemmamod_qbic_all %>% 
  filter(SNP %in% gemmamod_locus_snps_g_caQTL$SNP) %>% 
  filter(grepl("HES1", model))
