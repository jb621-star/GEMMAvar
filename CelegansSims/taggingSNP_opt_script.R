#library(genetics)
library(data.table)
library(vcfR)
library(tidyverse)
library(parallel)

setwd("/hpc/home/jb621")
# For the full data, see directory /hpc/group/baughlab/jb621/power_sims_scripts/EvoWorm2024_scripts
Garter_chrI_SNPs <- read.vcfR("/hpc/group/baughlab/jb621/power_sims_scripts/EvoWorm2024_scripts/WI.20220216.GarterSnake.I.SNPs_only.vcf.gz")
Garter_chrI_SNPs_df <- vcfR2tidy(Garter_chrI_SNPs, single_frame = TRUE, info_types = TRUE, format_types = TRUE)
Garter_chrI_SNPs_df_dat <- Garter_chrI_SNPs_df$dat
Garter_chrI_SNPs_df_Selalleles <- Garter_chrI_SNPs_df_dat %>%
  dplyr::select(CHROM, POS, REF, ALT, Indiv, gt_GT_alleles) %>%
  mutate(SNP = paste0(CHROM, ":", POS))

save(Garter_chrI_SNPs_df_Selalleles, file = "/work/jb621/Garter_chrI_SNPs_df_Selalleles.Rda")

load(file = "/work/jb621/simHaplo/LD_mapping/Garter_chrI_SNPs_df_Selalleles.Rda")

#load(file = "/work/jb621/simHaplo/LD_mapping/Garter_chrIV_SNPs_df_Selalleles.Rda")
#HapRepSNPs_chrI <- fread(file = "HapRepSNP_chr_I_data.tsv", sep = "\t")
#HapRepSNPs_chrI_gt <- HapRepSNPs_chrI %>% 
#  dplyr::select(V1,V2,V3) %>% 
#  pivot_wider(names_from = V2, values_from = V3)

Garter_chrI_SNPs_df_gt <- Garter_chrI_SNPs_df_Selalleles %>%
  dplyr::select(Indiv, SNP, gt_GT_alleles) %>%
  pivot_wider(names_from = SNP, values_from = gt_GT_alleles)

genotype_matrix <- as.matrix(Garter_chrI_SNPs_df_gt[,-1])
rownames(genotype_matrix) <- Garter_chrI_SNPs_df_gt$Indiv
#will give you the indices of columns
#cols <- genotype_matrix[, which(colMeans(genotype_matrix == ".") > 0.5)]
# cols <- c()
# for (i in 1:ncol(genotype_matrix)) {
#   if (sum(genotype_matrix[,i] == ".")/191 > 0.1) {
#     cols <- c(cols,i)
#   }
# }
# 
# genotype_matrix_missfilt <- genotype_matrix[,-cols]
# 
# sing <- c()
# for (i in 1:ncol(genotype_matrix)) {
#   if (length(unique(genotype_matrix_missfilt[,i])) != 3) {
#     sing <- c(sing,i)
#   }
# }
# Using vcfR to extract numeric genotype matrix
genotype_matrix_GF_SNP <- extract.gt(Garter_chrI_SNPs, as.numeric = T)
genotype_matrix_GF_SNP <- genotype_matrix_GF_SNP[complete.cases(genotype_matrix_GF_SNP), ]
genotype_matrix_GF_SNP <- t(genotype_matrix_GF_SNP)
# Change column names
colnames(genotype_matrix_GF_SNP) <- gsub("I_(\\d+)", "I:\\1", colnames(genotype_matrix_GF_SNP))

sing <- c()
for (i in 1:ncol(genotype_matrix_GF_SNP)) {
  if (length(unique(genotype_matrix_GF_SNP[,i])) != 2) {
    sing <- c(sing,i)
  }
}

genotype_matrix_GF_SNP_missfilt_bi <- genotype_matrix_GF_SNP[,-sing]
save(genotype_matrix_GF_SNP_missfilt_bi, file = "genotype_matrix_GF_SNP_missfilt_bi_chrI.Rda")
load(file = "genotype_matrix_GF_SNP_missfilt_bi_chrI.Rda")

snp_names <- colnames(genotype_matrix_missfilt_bi_chrI)

colnames(genotype_matrix_GF_SNP_missfilt_bi_chrI) <- gsub("I_(\\d+)", "I:\\1", colnames(genotype_matrix_GF_SNP_missfilt_bi_chrI))

genotype_matrix_GF_SNP_missfilt_bi_chrI_pt1 <- genotype_matrix_GF_SNP_missfilt_bi_chrI[,1:44724]
genotype_matrix_GF_SNP_missfilt_bi_chrI_pt2 <- genotype_matrix_GF_SNP_missfilt_bi_chrI[,44725:89447]

# genotype_matrix_GF_SNP_missfilt_bi_chrII_pt1 <- genotype_matrix_GF_SNP_missfilt_bi_chrII[,1:48328]
# genotype_matrix_GF_SNP_missfilt_bi_chrII_pt2 <- genotype_matrix_GF_SNP_missfilt_bi_chrII[,48329:96655]

# genotype_matrix_GF_SNP_missfilt_bi_chrIII_pt1 <- genotype_matrix_GF_SNP_missfilt_bi_chrIII[,1:39299]
# genotype_matrix_GF_SNP_missfilt_bi_chrIII_pt2 <- genotype_matrix_GF_SNP_missfilt_bi_chrIII[,39300:78598]

# genotype_matrix_GF_SNP_missfilt_bi_chrIV_pt1 <- genotype_matrix_GF_SNP_missfilt_bi_chrIV[,1:60882]
# genotype_matrix_GF_SNP_missfilt_bi_chrIV_pt2 <- genotype_matrix_GF_SNP_missfilt_bi_chrIV[,60883:121763]

# genotype_matrix_GF_SNP_missfilt_bi_chrX_pt1 <- genotype_matrix_GF_SNP_missfilt_bi_chrX[,1:69581]
# genotype_matrix_GF_SNP_missfilt_bi_chrX_pt2 <- genotype_matrix_GF_SNP_missfilt_bi_chrX[,69582:139161]

# genotype_matrix_GF_SNP_missfilt_bi_chrV_pt1 <- genotype_matrix_GF_SNP_missfilt_bi_chrV[,1:54898]
# genotype_matrix_GF_SNP_missfilt_bi_chrV_pt2 <- genotype_matrix_GF_SNP_missfilt_bi_chrV[,54899:109795]
# genotype_matrix_GF_SNP_missfilt_bi_chrV_pt3 <- genotype_matrix_GF_SNP_missfilt_bi_chrV[,109796:164693]
# genotype_matrix_GF_SNP_missfilt_bi_chrV_pt4 <- genotype_matrix_GF_SNP_missfilt_bi_chrV[,164694:219590]

# chrI_r2_mat_pt1 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrI_pt1)
# #save(chrI_r2_mat_pt1, file = "/work/jb621/simHaplo/LD_mapping/chrI_r2_mat_pt1.Rda")
# chrI_r2_mat_pt2 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrI_pt2)
chrI_r2_mat_pt1 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrI_pt1))^2
#save(chrII_r2_mat_pt1, file = "/work/jb621/simHaplo/LD_mapping/chrII_r2_mat_pt1.Rda")
#chrII_r2_mat_pt2 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrII_pt2)
chrI_r2_mat_pt2 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrI_pt2))^2

#chrII_r2_mat_pt1 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrII_pt1)
#chrII_r2_mat_pt1 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrII_pt1))^2
#save(chrII_r2_mat_pt1, file = "/work/jb621/simHaplo/LD_mapping/chrII_r2_mat_pt1.Rda")
#chrII_r2_mat_pt2 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrII_pt2)
chrII_r2_mat_pt2 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrII_pt2))^2

# chrIII_r2_mat_pt1 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrIII_pt1)
# #save(chrIII_r2_mat_pt1, file = "/work/jb621/simHaplo/LD_mapping/chrIII_r2_mat_pt1.Rda")
# chrIII_r2_mat_pt2 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrIII_pt2)
chrIII_r2_mat_pt1 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrIII_pt1))^2
#save(chrII_r2_mat_pt1, file = "/work/jb621/simHaplo/LD_mapping/chrII_r2_mat_pt1.Rda")
#chrII_r2_mat_pt2 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrII_pt2)
chrIII_r2_mat_pt2 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrIII_pt2))^2

chrIV_r2_mat_pt1 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrIV_pt1))^2
#save(chrII_r2_mat_pt1, file = "/work/jb621/simHaplo/LD_mapping/chrII_r2_mat_pt1.Rda")
#chrII_r2_mat_pt2 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrII_pt2)
chrIV_r2_mat_pt2 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrIV_pt2))^2

#chrV_r2_mat_pt1 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt1)
chrV_r2_mat_pt1 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt1))^2
#save(chrV_r2_mat_pt1, file = "/work/jb621/simHaplo/LD_mapping/chrV_r2_mat_pt1.Rda")
#chrV_r2_mat_pt2 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt2)
chrV_r2_mat_pt2 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt2))^2
#save(chrIV_r2_mat_pt2, file = "/work/jb621/simHaplo/LD_mapping/chrIV_r2_mat_pt2.Rda")
#chrV_r2_mat_pt3 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt3)
chrV_r2_mat_pt3 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt3))^2

#chrV_r2_mat_pt4 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt4)
chrV_r2_mat_pt4 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt4))^2


#chrX_r2_mat_pt1 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrX_pt1)
chrX_r2_mat_pt1 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrX_pt1))^2
#save(chrX_r2_mat_pt1, file = "/work/jb621/simHaplo/LD_mapping/chrX_r2_mat_pt1.Rda")
#chrX_r2_mat_pt2 <- calculate_r2_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrX_pt2)
chrX_r2_mat_pt2 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrX_pt2))^2

identify_haplotype_blocks_optimized <- function(genotype_matrix, r2_mat) {
  num_snps <- ncol(genotype_matrix)
  r2_matrix <- as.matrix(r2_mat)
  blocks <- list()
  used_snps <- rep(FALSE, num_snps)  # Boolean vector for quick lookups
  i <- 1
  
  while (i <= num_snps) {
    if (used_snps[i]) {
      i <- i + 1
      next  # Skip SNPs that are already used
    }
    
    current_block <- c(i)
    used_snps[i] <- TRUE
    j <- i + 1
    
    while (j <= num_snps) {
      if (used_snps[j]) {
        j <- j + 1
        next  # Skip SNPs that are already used
      }
      
      if (r2_matrix[i, j] >= 0.99) {
        current_block <- c(current_block, j)
        used_snps[j] <- TRUE
      }
      j <- j + 1
    }
    blocks[[length(blocks) + 1]] <- current_block
    i <- i + 1
  }
  
  return(blocks)
}
load(file = "genotype_matrix_missfilt_bi_II.Rda")

# genotype_matrix_missfilt_bi_pt1 <- genotype_matrix_missfilt_bi[,1:44724]
# genotype_matrix_missfilt_bi_pt2 <- genotype_matrix_missfilt_bi[,44725:89447]

# genotype_matrix_missfilt_bi_pt1 <- genotype_matrix_missfilt_bi[,1:48328]
# genotype_matrix_missfilt_bi_pt2 <- genotype_matrix_missfilt_bi[,48329:96655]

# genotype_matrix_missfilt_bi_pt1 <- genotype_matrix_missfilt_bi[,1:39299]
# genotype_matrix_missfilt_bi_pt2 <- genotype_matrix_missfilt_bi[,39300:78598]

# genotype_matrix_missfilt_bi_pt1 <- genotype_matrix_missfilt_bi[,1:60882]
# genotype_matrix_missfilt_bi_pt2 <- genotype_matrix_missfilt_bi[,60883:121763]

# genotype_matrix_missfilt_bi_pt1 <- genotype_matrix_missfilt_bi[,1:54898]
# genotype_matrix_missfilt_bi_pt2 <- genotype_matrix_missfilt_bi[,54899:109795]
# genotype_matrix_missfilt_bi_pt3 <- genotype_matrix_missfilt_bi[,109796:164693]
# genotype_matrix_missfilt_bi_pt4 <- genotype_matrix_missfilt_bi[,164694:219590]

# genotype_matrix_missfilt_bi_pt1 <- genotype_matrix_missfilt_bi[,1:69581]
# genotype_matrix_missfilt_bi_pt2 <- genotype_matrix_missfilt_bi[,69582:139161]

snp_names_pt1 <- colnames(genotype_matrix_GF_SNP_missfilt_bi_chrX_pt1)
snp_names_pt2 <- colnames(genotype_matrix_GF_SNP_missfilt_bi_chrX_pt2)
# snp_names_pt3 <- colnames(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt3)
# snp_names_pt4 <- colnames(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt4)

# Function to convert snp_groups from indices to SNP names
convert_snp_groups_to_names <- function(snp_groups, snp_names) {
  # Initialize an empty list to store the converted groups
  snp_groups_named <- lapply(snp_groups, function(group) {
    # Map the indices in each group to the corresponding SNP names
    snp_names[group]
  })
  
  return(snp_groups_named)
}

chrI_pt1_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrI_pt1, chrI_r2_mat_pt1_v2)
save(chrI_pt1_blocks_v2, file = "chrI_withlowMAF_pt1_blocks_v2.Rda")
load(file = "chrI_withlowMAF_pt1_blocks.Rda")

sum(lengths(chrI_pt1_blocks))

chrI_pt1_blocks_named_v2 <- convert_snp_groups_to_names(chrI_pt1_blocks_v2, snp_names_pt1)

chrI_pt2_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrI_pt2, chrI_r2_mat_pt2_v2)
save(chrI_pt2_blocks_v2, file = "chrI_withlowMAF_pt2_blocks_v2.Rda")
load(file = "chrI_withlowMAF_pt2_blocks.Rda")

sum(lengths(chrI_pt2_blocks))

chrI_pt2_blocks_named_v2 <- convert_snp_groups_to_names(chrI_pt2_blocks_v2, snp_names_pt2)

chrI_all_blocks_named_v2 <- c(chrI_pt1_blocks_named_v2, chrI_pt2_blocks_named_v2)
save(chrI_all_blocks_named_v2, file = "chrI_withlowMAF_all_blocks_named_v2.Rda")

#chrII

chrII_pt1_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrII_pt1, chrII_r2_mat_pt1_v2)
save(chrII_pt1_blocks_v2, file = "chrII_withlowMAF_pt1_blocks_v2.Rda")
load(file = "chrII_withlowMAF_pt1_blocks_v2.Rda")

chrII_pt1_blocks_named_v2 <- convert_snp_groups_to_names(chrII_pt1_blocks_v2, snp_names_pt1)

chrII_pt1_tagging_SNPs_v2 <- sapply(chrII_pt1_blocks_named_v2, `[[`, 1)


# Find specific SNPs...why are some perfect LD SNPs still in here?
# Find index of element containing "II:6081257"
# idx <- which(sapply(chrII_pt1_blocks_named, function(x) "II:6467210" %in% x))
# 
# # Extract the matching element
# matching_element <- chrII_pt1_blocks_named[[idx]]

chrII_pt2_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrII_pt2, chrII_r2_mat_pt2_v2)
save(chrII_pt2_blocks_v2, file = "chrII_withlowMAF_pt2_blocks_v2.Rda")
load(file = "chrII_withlowMAF_pt2_blocks_v2.Rda")

chrII_pt2_blocks_named_v2 <- convert_snp_groups_to_names(chrII_pt2_blocks_v2, snp_names_pt2)

chrII_pt2_tagging_SNPs_v2 <- sapply(chrII_pt2_blocks_named_v2, `[[`, 1)

chrII_pretagging_SNPs_v2 <- c(chrII_pt1_tagging_SNPs_v2, chrII_pt2_tagging_SNPs_v2)

genotype_matrix_GF_SNP_missfilt_bi_chrII_pretagging <- genotype_matrix_GF_SNP_missfilt_bi_chrII[,chrII_pretagging_SNPs_v2]

chrII_r2_mat_pretagging_v2 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrII_pretagging))^2

chrII_all_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrII_pretagging, chrII_r2_mat_pretagging_v2)
save(chrII_all_blocks_v2, file = "chrII_withlowMAF_all_blocks_v2.Rda")

# snp_names_pretagging <- colnames(genotype_matrix_GF_SNP_missfilt_bi_chrII_pretagging)
# 
# chrII_all_blocks_named_v2 <- convert_snp_groups_to_names(chrII_all_blocks_v2, snp_names_pretagging)

chrII_all_blocks_named_v2 <- c(chrII_pt1_blocks_named_v2, chrII_pt2_blocks_named_v2)
save(chrII_all_blocks_named_v2, file = "chrII_withlowMAF_all_blocks_named_v2.Rda")

# Find specific SNPs...why are some perfect LD SNPs still in here?
# Find index of element containing "II:6081257"
# idx <- which(sapply(chrII_pt2_blocks_named_v2, function(x) "II:10049635" %in% x))
# 
# # Extract the matching element
# matching_element <- chrII_pt2_blocks_named_v2[[idx]]

# chrII_all_blocks_named <- c(chrII_pt1_blocks_named, chrII_pt2_blocks_named)
# save(chrII_all_blocks_named, file = "chrII_withlowMAF_all_blocks_named.Rda")

# no missing entries...but somehow 

#chrIII

chrIII_pt1_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrIII_pt1, chrIII_r2_mat_pt1_v2)
save(chrIII_pt1_blocks_v2, file = "chrIII_withlowMAF_pt1_blocks_v2.Rda")
#load(file = "chrIII_withlowMAF_pt1_blocks.Rda")

chrIII_pt1_blocks_named_v2 <- convert_snp_groups_to_names(chrIII_pt1_blocks_v2, snp_names_pt1)

chrIII_pt2_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrIII_pt2, chrIII_r2_mat_pt2_v2)
save(chrIII_pt2_blocks_v2, file = "chrIII_withlowMAF_pt2_blocks_v2.Rda")
#load(file = "chrIII_withlowMAF_pt2_blocks.Rda")

chrIII_pt2_blocks_named_v2 <- convert_snp_groups_to_names(chrIII_pt2_blocks_v2, snp_names_pt2)

chrIII_all_blocks_named_v2 <- c(chrIII_pt1_blocks_named_v2, chrIII_pt2_blocks_named_v2)
save(chrIII_all_blocks_named_v2, file = "chrIII_withlowMAF_all_blocks_named_v2.Rda")

# chrIII_all_blocks <- c(chrIII_pt1_blocks, chrIII_pt2_blocks)
# save(chrIII_all_blocks, file = "chrIII_withlowMAF_all_blocks.Rda")

# chrIV

chrIV_pt1_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrIV_pt1, chrIV_r2_mat_pt1_v2)
save(chrIV_pt1_blocks_v2, file = "chrIV_withlowMAF_pt1_blocks_v2.Rda")
#load(file = "chrIV_withlowMAF_pt1_blocks.Rda")

chrIV_pt1_blocks_named_v2 <- convert_snp_groups_to_names(chrIV_pt1_blocks_v2, snp_names_pt1)

chrIV_pt2_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrIV_pt2, chrIV_r2_mat_pt2_v2)
save(chrIV_pt2_blocks_v2, file = "chrIV_withlowMAF_pt2_blocks_v2.Rda")
#load(file = "chrIV_withlowMAF_pt2_blocks.Rda")

chrIV_pt2_blocks_named_v2 <- convert_snp_groups_to_names(chrIV_pt2_blocks_v2, snp_names_pt2)

chrIV_all_blocks_named_v2 <- c(chrIV_pt1_blocks_named_v2, chrIV_pt2_blocks_named_v2)
save(chrIV_all_blocks_named_v2, file = "chrIV_withlowMAF_all_blocks_named_v2.Rda")

# chrIV_all_blocks <- c(chrIV_pt1_blocks, chrIV_pt2_blocks)
# save(chrIV_all_blocks, file = "chrIV_withlowMAF_all_blocks.Rda")

# chrV
chrV_pt1_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt1, chrV_r2_mat_pt1_v2)
save(chrV_pt1_blocks_v2, file = "chrV_withlowMAF_pt1_blocks_v2.Rda")
#load(file = "chrV_withlowMAF_pt1_blocks.Rda")

chrV_pt1_blocks_named_v2 <- convert_snp_groups_to_names(chrV_pt1_blocks_v2, snp_names_pt1)

chrV_pt2_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt2, chrV_r2_mat_pt2_v2)
save(chrV_pt2_blocks_v2, file = "chrV_withlowMAF_pt2_blocks_v2.Rda")
#load(file = "chrV_withlowMAF_pt2_blocks.Rda")

chrV_pt2_blocks_named_v2 <- convert_snp_groups_to_names(chrV_pt2_blocks_v2, snp_names_pt2)

chrV_pt3_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt3, chrV_r2_mat_pt3_v2)
save(chrV_pt3_blocks_v2, file = "chrV_withlowMAF_pt3_blocks_v2.Rda")
#load(file = "chrV_withlowMAF_pt3_blocks.Rda")

chrV_pt3_blocks_named_v2 <- convert_snp_groups_to_names(chrV_pt3_blocks_v2, snp_names_pt3)

chrV_pt4_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrV_pt4, chrV_r2_mat_pt4_v2)
save(chrV_pt4_blocks_v2, file = "chrV_withlowMAF_pt4_blocks_v2.Rda")
#load(file = "chrV_withlowMAF_pt4_blocks.Rda")

chrV_pt4_blocks_named_v2 <- convert_snp_groups_to_names(chrV_pt4_blocks_v2, snp_names_pt4)

chrV_all_blocks_named_v2 <- c(chrV_pt1_blocks_named_v2, chrV_pt2_blocks_named_v2, chrV_pt3_blocks_named_v2, chrV_pt4_blocks_named_v2)
save(chrV_all_blocks_named_v2, file = "chrV_withlowMAF_all_blocks_named_v2.Rda")
#load(file = "chrV_withlowMAF_all_blocks_named_v2.Rda")

# chrX
# At once version
chrX_r2_mat_atonce_v2 <- cor(genotype_matrix_GF_SNP_missfilt_bi_chrX)

chrX_pt1_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrX_pt1, chrX_r2_mat_pt1_v2)
save(chrX_pt1_blocks_v2, file = "chrX_withlowMAF_pt1_blocks_v2.Rda")
#load(file = "chrX_withlowMAF_pt1_blocks.Rda")

chrX_pt1_blocks_named_v2 <- convert_snp_groups_to_names(chrX_pt1_blocks_v2, snp_names_pt1)

chrX_pt2_blocks_v2 <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrX_pt2, chrX_r2_mat_pt2_v2)
save(chrX_pt2_blocks_v2, file = "chrX_withlowMAF_pt2_blocks_v2.Rda")
#load(file = "chrX_withlowMAF_pt2_blocks.Rda")

chrX_pt2_blocks_named_v2 <- convert_snp_groups_to_names(chrX_pt2_blocks_v2, snp_names_pt2)

chrX_all_blocks_named_v2 <- c(chrX_pt1_blocks_named_v2, chrX_pt2_blocks_named_v2)
chrX_all_blocks_named_v2 <- lapply(chrX_all_blocks_named_v2, function(x) gsub("_", ":", x))
save(chrX_all_blocks_named_v2, file = "chrX_withlowMAF_all_blocks_named_v2.Rda")

# place all blocks together and save
# load(file = "chrI_withlowMAF_all_blocks_named_v2.Rda")
# load(file = "chrII_withlowMAF_all_blocks_named_v2.Rda")
# load(file = "chrIII_withlowMAF_all_blocks_named_v2.Rda")
#load(file = "chrIV_withlowMAF_all_blocks_named_v2.Rda")
#load(file = "chrV_withlowMAF_all_blocks_named_v2.Rda")
# load(file = "chrI_withlowMAF_all_blocks_named_v2.Rda")

load(file = "chrI_withlowMAF_atonce_blocks_named_v2.Rda")
load(file = "chrII_withlowMAF_atonce_blocks_named_v2.Rda")
load(file = "chrIII_withlowMAF_atonce_blocks_named_v2.Rda")
load(file = "chrIV_withlowMAF_atonce_blocks_named_v2.Rda")
load(file = "chrV_withlowMAF_atonce_blocks_named_v2.Rda")
load(file = "chrX_withlowMAF_atonce_blocks_named_v2.Rda")

allchr_blocks_named_v2 <- c(chrI_blocks_named_v2, chrII_blocks_named_v2,
                            chrIII_blocks_named_v2, chrIV_blocks_named_v2,
                            chrV_blocks_named_v2, chrX_blocks_named_v2)

# allchr_blocks_named_v2 <- c(chrI_all_blocks_named_v2, chrII_all_blocks_named_v2,
#                             chrIII_all_blocks_named_v2, chrIV_all_blocks_named_v2,
#                             chrV_all_blocks_named_v2, chrX_all_blocks_named_v2)

#save(allchr_blocks_named_v2, file = "allchr_blocks_atonce_named_v2.Rda")
load(file = "allchr_blocks_atonce_named_v2.Rda")
# This was to make sure that there are no singleton tagging SNPs, v2 --> v3
# Filter the list based on the first element of each vector
allchr_blocks_named_v3 <- allchr_blocks_named_v2[sapply(allchr_blocks_named_v2, function(x) x[1] %in% all_tagging_SNPs)]

save(allchr_blocks_named_v3, file = "allchr_blocks_atonce_named_v3.Rda")

prune_snps_in_blocks <- function(blocks, snp_names) {
  #start_time <- Sys.time()
  rep_snps_all <- c()
  for (i in 1:length(blocks)) {
    block_snps <- blocks[[i]]
    rep_snp <- block_snps[1] # Keep the first SNP in each block
    rep_snp_name <- snp_names[rep_snp]
    rep_snps_all[i] <- rep_snp_name
    
    #for (snp in block[-1]) {
    #  if ((snp_positions[snp] - last_position) >= (1 / cM_per_Mb) * 1e6 * 50) {
    #    pruned_snps <- c(pruned_snps, snp)
    #    last_position <- snp_positions[snp]
    #  }
    #}
  }
  #end_time <- Sys.time()
  #print(paste("prune_snps_in_blocks function took", end_time - start_time, "seconds"))
  
  return(rep_snps_all)
}

chrI_HapRepSNPs_pt1 <- prune_snps_in_blocks(chrI_pt1_blocks, snp_names = snp_names_pt1)
chrI_HapRepSNPs_pt2 <- prune_snps_in_blocks(chrI_pt2_blocks, snp_names = snp_names_pt2)
chrI_HapRepSNPs_rnd2 <- c(chrI_HapRepSNPs_pt1, chrI_HapRepSNPs_pt2)

colnames(genotype_matrix_GF_SNP_missfilt_bi_chrI) <- gsub("I_(\\d+)", "I:\\1", colnames(genotype_matrix_GF_SNP_missfilt_bi_chrI))
genotype_matrix_GF_SNP_missfilt_bi_chrI_rnd2 <- genotype_matrix_GF_SNP_missfilt_bi_chrI[,colnames(genotype_matrix_GF_SNP_missfilt_bi_chrI) %in% chrI_HapRepSNPs_rnd2]

chrI_r2_mat_rnd2 <- (cor(genotype_matrix_GF_SNP_missfilt_bi_chrI_rnd2)^2)


chrI_rnd2_blocks <- identify_haplotype_blocks_optimized(genotype_matrix_GF_SNP_missfilt_bi_chrI_rnd2, chrI_r2_mat_rnd2)
save(chrI_rnd2_blocks, file = "chrI_withlowMAF_rnd2_blocks.Rda")
load(file = "chrI_withlowMAF_pt1_blocks.Rda")
load(file = "chrI_withlowMAF_pt2_blocks.Rda")
chrI_all_blocks <- c(chrI_pt1_blocks, chrI_pt2_blocks)
save(chrI_all_blocks, file = "chrI_withlowMAF_all_blocks.Rda")

load(file = "chrI_withlowMAF_all_blocks.Rda")

calculate_maf <- function(genotype_column) {
  allele_counts <- table(genotype_column)
  total_alleles <- sum(allele_counts)
  
  # If the genotype is coded as 0, 1, 2 for number of alternative alleles
  min_allele_count <- min(allele_counts[1], allele_counts[2], na.rm = TRUE)
  maf <- min_allele_count / total_alleles
  if (maf == 1) {
    maf <- 0
  }
  return(maf)
}

# Function to filter SNPs based on MAF
filter_snps_by_maf <- function(genotype_matrix, maf_threshold = 0.01) {
  maf_values <- apply(genotype_matrix, 2, calculate_maf)
  filtered_genotype_matrix <- genotype_matrix[, maf_values > maf_threshold]
  return(filtered_genotype_matrix)
}

MAFfiltered_genotypes_chrI <- filter_snps_by_maf(genotype_matrix_missfilt_bi)
#save(MAFfiltered_genotypes_chrI, file = "genotype_matrix_missfilt_bi_maf_I.Rda")


#load(file = "genotype_matrix_missfilt_bi_maf_I.Rda")
#load(file = "genotype_matrix_missfilt_bi_maf_X.Rda")
load(file = "genotype_matrix_missfilt_bi_I.Rda")

