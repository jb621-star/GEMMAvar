# Load necessary library
library(data.table)
library(susieR)
library(tidyverse)
#library(RcppGSL)
library(Rfast)
library(curl)
library(vcfR)

load(file = "allchr_blocks_atonce_named_v3.Rda")

#save(NORMAL.genomat.filt_mat, file = "/work/jb621/simHaplo/LD_mapping/taggingSNPs_v4_geno_df.Rda")
load(file = "NORMgenomat.samp.Rda")
load(file = "taggingSNPs_v4_geno_df.Rda")

taggingSNPs_v4_geno_df <- NORMAL.genomat.filt_mat / 2 

library(dplyr)

split_snps <- strsplit(colnames(taggingSNPs_v4_geno_df), ":")

# Create a dataframe by applying a function to split_snps
snp_info <- do.call(rbind, lapply(split_snps, function(snp) {
  chromosome <- snp[1]                 # Chromosome code (e.g., "I")
  start_range <- as.numeric(snp[2])     # Start of range (numeric position)
  
  set_id <- paste0(chromosome, ":", start_range)  # Set ID
  
  # Return as a row in the dataframe
  return(data.frame(chromosome, set_id))
}))

summarize_snp_maf_ld <- function(genotype_matrix, snp_info) {
  genotype_matrix <- taggingSNPs_v4_geno_df
  
  # Step 1: Calculate MAF
  allele_freqs <- colMeans(genotype_matrix, na.rm = T)  # assuming 0/1 coding
  mafs <- pmin(allele_freqs, 1 - allele_freqs)
  
  # Step 2: Standardize genotypes
  centered_genotype_matrix <- scale(genotype_matrix, center = TRUE, scale = TRUE)
  
  # Step 3: Compute LD matrix (R²)
  ld_matrix <- cor(centered_genotype_matrix, use = "pairwise.complete.obs")^2
  
  #genotype_matrix <- t(genotype_matrix)
  n_snps <- ncol(genotype_matrix)
  snp_ids <- colnames(genotype_matrix)
  
  # Step 4: Initialize result vectors
  avg_cis_ld <- numeric(n_snps)
  med_cis_ld <- numeric(n_snps)
  max_cis_ld <- numeric(n_snps)
  avg_trans_ld <- numeric(n_snps)
  med_trans_ld <- numeric(n_snps)
  max_trans_ld <- numeric(n_snps)
  
  for (i in 1:n_snps) {
    #i <- 1
    chr_i <- snp_info$chromosome[i]
    
    cis_indices <- which(snp_info$chromosome == chr_i)
    trans_indices <- which(snp_info$chromosome != chr_i)
    cis_indices <- setdiff(cis_indices, i)  # exclude self
    
    cis_ld_vals <- if (length(cis_indices) > 0) ld_matrix[i, cis_indices] else NA
    trans_ld_vals <- if (length(trans_indices) > 0) ld_matrix[i, trans_indices] else NA
    
    avg_cis_ld[i]   <- mean(cis_ld_vals, na.rm = TRUE)
    med_cis_ld[i]   <- median(cis_ld_vals, na.rm = TRUE)
    max_cis_ld[i]   <- max(cis_ld_vals, na.rm = TRUE)
    avg_trans_ld[i] <- mean(trans_ld_vals, na.rm = TRUE)
    med_trans_ld[i] <- median(trans_ld_vals, na.rm = TRUE)
    max_trans_ld[i] <- max(trans_ld_vals, na.rm = TRUE)
  }
  
  # Step 5: Combine into a data frame
  snp_summary <- data.frame(
    SNP = snp_ids,
    Chromosome = snp_info$chromosome,
    MAF = mafs,
    Avg_Cis_LD = avg_cis_ld,
    Med_Cis_LD = med_cis_ld,
    Max_Cis_LD = max_cis_ld,
    Avg_Trans_LD = avg_trans_ld,
    Med_Trans_LD = med_trans_ld,
    Max_Trans_LD = max_trans_ld,
    stringsAsFactors = FALSE
  )
  
  return(snp_summary)
}

#taggingSNPs_v3_summary <- summarize_snp_maf_ld(taggingSNPs_v3_geno_df, snp_info)
taggingSNPs_v4_summary <- snp_summary
save(taggingSNPs_v4_summary, file = "/work/jb621/simHaplo/LD_mapping/taggingSNPs_v4_summary.Rda")
#load(file = "/work/jb621/simHaplo/LD_mapping/taggingSNPs_v4_summary.Rda")

# Random sample of 1000 SNPs, then analyze the MAF, LD context of them
set.seed(43)
taggingSNPs_v4_geno_df <- data.frame(t(taggingSNPs_v4_geno_df))
dim(taggingSNPs_v4_geno_df)
sampled.taggingSNPs <- taggingSNPs_v4_geno_df %>% 
  slice_sample(n=1000)

taggingSNPs_v4_summary_sampled <- taggingSNPs_v4_summary %>% 
  filter(SNP %in% rownames(sampled.taggingSNPs))

library(ggExtra)

# Plot relationship between avg and max
r <- cor(taggingSNPs_v4_summary$Avg_Cis_LD, taggingSNPs_v4_summary$Max_Cis_LD, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))
IntraScatter <- ggplot(taggingSNPs_v4_summary, aes(x = Avg_Cis_LD, y = Max_Cis_LD)) +
  geom_point() +
  labs(x = "Avg Intra Coinheritance", y = "Max Intra Coinheritance", title = "Relationship between Intra Coinheritance Measures") +
  theme_minimal() + 
  annotate("text", x = Inf, y = 0.25, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18))             # Increase y-axis label size
ggMarginal(IntraScatter, type = "histogram")

r <- cor(taggingSNPs_v4_summary$Avg_Cis_LD, taggingSNPs_v4_summary$Med_Cis_LD, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))
IntraCenterScatter <- ggplot(taggingSNPs_v4_summary, aes(x = Avg_Cis_LD, y = Med_Cis_LD)) +
  geom_point() +
  labs(x = "Avg Intra Coinheritance", y = "Median Intra Coinheritance", title = "Relationship between Intra Coinheritance Measures") +
  theme_minimal() + 
  annotate("text", x = 0.025, y = 0.04, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18)) 
ggMarginal(IntraCenterScatter, type = "histogram")

r <- cor(taggingSNPs_v4_summary$Avg_Trans_LD, taggingSNPs_v4_summary$Max_Trans_LD, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))
InterScatter <- ggplot(taggingSNPs_v4_summary, aes(x = Avg_Trans_LD, y = Max_Trans_LD)) +
  geom_point() +
  labs(x = "Avg Inter Coinheritance", y = "Max Inter Coinheritance", title = "Relationship between Inter Coinheritance Measures") +
  theme_minimal() + 
  annotate("text", x = 0.045, y = 0.375, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18))             # Increase y-axis label size
ggMarginal(InterScatter, type = "histogram")

r <- cor(taggingSNPs_v4_summary$Avg_Trans_LD, taggingSNPs_v4_summary$Med_Trans_LD, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))
InterCenterScatter <- ggplot(taggingSNPs_v4_summary, aes(x = Avg_Trans_LD, y = Med_Trans_LD)) +
  geom_point() +
  labs(x = "Avg Inter Coinheritance", y = "Median Inter Coinheritance", title = "Relationship between Inter Coinheritance Measures") +
  theme_minimal() + 
  annotate("text", x = 0.05, y = 0.0125, label = r2_label,
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18)) 
ggMarginal(InterCenterScatter, type = "histogram")

taggingSNPs_v4_summary %>% 
  reframe(quantile = scales::percent(c(0.25, 0.5, 0.75)),
          MAF = quantile(MAF, c(0.25, 0.5, 0.75)),
          Avg_Cis_LD = quantile(Avg_Cis_LD, c(0.25, 0.5, 0.75)),
          Med_Cis_LD = quantile(Med_Cis_LD, c(0.25, 0.5, 0.75)),
          Max_Cis_LD = quantile(Max_Cis_LD, c(0.25, 0.5, 0.75)),
          Avg_Trans_LD = quantile(Avg_Trans_LD, c(0.25, 0.5, 0.75)),
          Med_Trans_LD = quantile(Med_Trans_LD, c(0.25, 0.5, 0.75)),
          Max_Trans_LD = quantile(Max_Trans_LD, c(0.25, 0.5, 0.75)))

# Examining the distribution of these characteristics in the sample, I chose the following
# as thresholds: 
taggingSNPs_v4_summary_sampled_cat <- taggingSNPs_v4_summary_sampled %>% 
  summarise(n_lowMAF = sum(MAF < 0.25),
            n_highMAF = sum(MAF >= 0.25),
            n_lowAvgIntraCoin = sum(Avg_Cis_LD < 0.05),
            n_lowMaxIntraCoin = sum(Max_Cis_LD < 0.5),
            n_highAvgIntraCoin = sum(Avg_Cis_LD > 0.05),
            n_highMaxIntraCoin = sum(Max_Cis_LD > 0.5),
            n_lowAvgInterCoin = sum(Avg_Trans_LD < 0.025),
            n_lowMaxInterCoin = sum(Max_Trans_LD < 0.5),
            n_highAvgInterCoin = sum(Avg_Trans_LD > 0.025),
            n_highMaxInterCoin = sum(Max_Trans_LD > 0.5))
# According to our categories, the sample has:
# majority low MAF SNPs

# majority low AVERAGE intra-coinheritance
# majority high MAX intra-coinheritance

# majority low AVERAGE inter-coinheritance
# majority low MAX inter-coinheritance

# Based on this initial survey, sampling needs to cover for gaps in:
# high MAF SNPs

# high AVERAGE intra-coinheritance
# low MAX intra-coinheritance

# high AVERAGE inter-coinheritance
# high MAX inter-coinheritance

# Inverse eCDF method 
# What exactly am I fitting the function to? MaxIntra, MaxInter, MAF separately?
# And then sampling roughly equally from those three distributions? 
# There does tend to be a correlation between R2 and MAF according to Paul, but is
# that what we observe?
r <- cor(taggingSNPs_v4_summary$Max_Cis_LD, taggingSNPs_v4_summary$MAF, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))

MaxIntraMAFScatter <- ggplot(taggingSNPs_v4_summary, aes(x = Max_Cis_LD, y = MAF)) +
  geom_point() +
  labs(
    x = "Max Intra Coinheritance", 
    y = "MAF", 
    title = "Relationship between Intra Coinheritance and MAF"
  ) +
  theme_minimal() + 
  annotate("text", x = 0.24, y = Inf, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

MaxIntraMAFScatter
ggMarginal(MaxIntraMAFScatter, type = "histogram")

r <- cor(taggingSNPs_v4_summary$Avg_Cis_LD, taggingSNPs_v4_summary$MAF, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))

AvgIntraMAFScatter <- ggplot(taggingSNPs_v4_summary, aes(x = Avg_Cis_LD, y = MAF)) +
  geom_point() +
  labs(
    x = "Avg Intra Coinheritance", 
    y = "MAF", 
    title = "Relationship between Avg Intra Coinheritance and MAF"
  ) +
  theme_minimal() + 
  annotate("text", x = 0.125, y = Inf, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

AvgIntraMAFScatter
ggMarginal(AvgIntraMAFScatter, type = "histogram")


r <- cor(taggingSNPs_v4_summary$Max_Trans_LD, taggingSNPs_v4_summary$MAF, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))

MaxInterMAFScatter <- ggplot(taggingSNPs_v4_summary, aes(x = Max_Trans_LD, y = MAF)) +
  geom_point() +
  labs(
    x = "Max Inter Coinheritance", 
    y = "MAF", 
    title = "Relationship between Inter Coinheritance and MAF"
  ) +
  theme_minimal() + 
  annotate("text", x = Inf, y = Inf, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

MaxInterMAFScatter
ggMarginal(MaxInterMAFScatter, type = "histogram")

r <- cor(taggingSNPs_v4_summary$Avg_Trans_LD, taggingSNPs_v4_summary$MAF, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))

AvgInterMAFScatter <- ggplot(taggingSNPs_v4_summary, aes(x = Avg_Trans_LD, y = MAF)) +
  geom_point() +
  labs(
    x = "Avg Inter Coinheritance", 
    y = "MAF", 
    title = "Relationship between Avg Inter Coinheritance and MAF"
  ) +
  theme_minimal() + 
  annotate("text", x = Inf, y = Inf, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

AvgInterMAFScatter
ggMarginal(AvgInterMAFScatter, type = "histogram")


r <- cor(taggingSNPs_v4_summary$Max_Trans_LD, taggingSNPs_v4_summary$Max_Cis_LD, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))
MaxCoinScatter <- ggplot(taggingSNPs_v4_summary, aes(x = Max_Trans_LD, y = Max_Cis_LD)) +
  geom_point() +
  labs(x = "Max Inter Coinheritance", y = "Max Intra Coinheritance", title = "Relationship between Inter and Intra Coinheritance") +
  theme_minimal() + 
  annotate("text", x = Inf, y = Inf, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18)) 
MaxCoinScatter
ggMarginal(MaxCoinScatter, type = "histogram")

# Need to derive a similar set of SNPs for the full 550 strains (at some point)
load(file = "/work/jb621/simHaplo/LD_mapping/taggingSNPs_v4_summary.Rda")

# Using the DENSITY function, not the ecdf function
# Compute empirical density estimate
density_est <- density(taggingSNPs_v4_summary$Max_Cis_LD, from = min(taggingSNPs_v4_summary$Max_Cis_LD), to = max(taggingSNPs_v4_summary$Max_Cis_LD))

# Interpolate density values at each data point
density_values <- approx(density_est$x, density_est$y, xout = taggingSNPs_v4_summary$Max_Cis_LD, rule = 2)$y

# Invert the density: higher weights for rare values
sampling_weights <- 1 / density_values

# Normalize weights
sampling_weights <- sampling_weights / sum(sampling_weights)

# Draw weighted sample
set.seed(43)
sample_indices <- sample(seq_along(sampling_weights), size = 1000, replace = FALSE, prob = sampling_weights)
x_samples_uniform <- taggingSNPs_v4_summary$Max_Cis_LD[sample_indices]

# Plot to confirm
p1 <- hist(taggingSNPs_v4_summary$Max_Cis_LD, breaks = 40, freq = FALSE)
p2 <- hist(x_samples_uniform, breaks = 40, xlab = "Distribution of Max Intra Coinheritance among sampled SNPs", main = "invPDF sampling based on Max Intra Coinheritance")
#plot( p1, col=rgb(0,0,1,1/4), xlab = "Distribution of Max Inter Coinheritance among tagging SNPs", main = "Original SNP distribution vs Weighted-Uniform Sample based on PDF")  # first histogram
plot( p1, col=rgb(0,0,1,1/4), xlab = "Distribution of Max Intra Coinheritance among tagging SNPs", main = "Original SNP distribution vs Weighted-Uniform Sample based on PDF")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), add=T)  # second
legend("topleft", legend = c("Original", "Uniform-sampled"), col = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), lwd = 2)

taggingSNPs_v4_summary_invpdf <- taggingSNPs_v4_summary[sample_indices,]

taggingSNPs_v4_summary_invpdf <- taggingSNPs_v4_summary_invpdf %>% 
  separate(SNP, into = c("chr", "pos"), sep = ":") %>% 
  dplyr::select(!Chromosome) %>% 
  mutate(pos = as.numeric(pos)) %>% 
  arrange(chr, pos) 
save(taggingSNPs_v4_summary_invpdf, file = "/work/jb621/simHaplo/LD_mapping/ars_SS_u7o3_sims/StepWiseSuSiE/invpdf_sampling_causal/taggingSNPs_v4_summary_invpdf.Rda")
load(file = "/work/jb621/simHaplo/LD_mapping/ars_SS_u7o3_sims_v2/ars_SS_u7o3_sims/StepWiseSuSiE/invpdf_sampling_causal/taggingSNPs_v4_summary_invpdf.Rda")

# Just as of curiosity, how well have we covered the whole genome?
taggingSNPs_v4_summary_invpdf %>% 
  group_by(chr) %>% 
  summarise(n = n()) # ended up being pretty evenly distributed throughout the genome
# selecting around 150 snps per chromosome

# Now checking to see if the distribution mimics the population along the other dimensions
# Compute R²
r <- cor(taggingSNPs_v4_summary_invpdf$Max_Cis_LD, taggingSNPs_v4_summary_invpdf$MAF, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))

library(ggExtra)
MaxIntraMAFScatter_sample <- ggplot(taggingSNPs_v4_summary_invpdf, aes(x = Max_Cis_LD, y = MAF)) +
  geom_point() +
  labs(x = "Max Intra Coinheritance", y = "MAF", title = "Relationship between Intra Coinheritance and MAF within sample") +
  theme_minimal() + 
  annotate("text", x = 0.3, y = Inf, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18)) 
MaxIntraMAFScatter_sample
ggMarginal(MaxIntraMAFScatter_sample, type = "histogram")
ggMarginal(MaxIntraMAFScatter, type = "histogram")

r <- cor(taggingSNPs_v4_summary_invpdf$Max_Trans_LD, taggingSNPs_v4_summary_invpdf$MAF, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))

MaxInterMAFScatter_sample <- ggplot(taggingSNPs_v4_summary_invpdf, aes(x = Max_Trans_LD, y = MAF)) +
  geom_point() +
  labs(x = "Max Inter Coinheritance", y = "MAF", title = "Relationship between Inter Coinheritance and MAF within sample") +
  theme_minimal() + 
  annotate("text", x = Inf, y = Inf, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18)) 
MaxInterMAFScatter_sample
ggMarginal(MaxInterMAFScatter_sample, type = "histogram")
ggMarginal(MaxInterMAFScatter, type = "histogram")

r <- cor(taggingSNPs_v4_summary_invpdf$Max_Trans_LD, taggingSNPs_v4_summary_invpdf$Max_Cis_LD, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))

MaxCoinScatter_sample <- ggplot(taggingSNPs_v4_summary_invpdf, aes(x = Max_Trans_LD, y = Max_Cis_LD)) +
  geom_point() +
  labs(x = "Max Inter Coinheritance", y = "Max Intra Coinheritance", title = "Relationship between Inter and Intra Coinheritance within sample") +
  theme_minimal() + 
  annotate("text", x = Inf, y = Inf, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18)) 
MaxCoinScatter_sample
ggMarginal(MaxCoinScatter_sample, type = "histogram")
ggMarginal(MaxCoinScatter, type = "histogram")

r <- cor(taggingSNPs_v4_summary_invpdf$Avg_Cis_LD, taggingSNPs_v4_summary_invpdf$Max_Cis_LD, use = "complete.obs")
r2 <- r^2
r2_label <- paste0("R2 = ", list(r2_val = format(r2, digits = 3)))
IntraScatter <- ggplot(taggingSNPs_v4_summary_invpdf, aes(x = Avg_Cis_LD, y = Max_Cis_LD)) +
  geom_point() +
  labs(x = "Avg Intra Coinheritance", y = "Max Intra Coinheritance", title = "Relationship between Intra Coinheritance Measures within sample") +
  theme_minimal() + 
  annotate("text", x = Inf, y = 0.25, label = r2_label, 
           hjust = 1.1, vjust = 1.5, size = 6) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18))             # Increase y-axis label size
ggMarginal(IntraScatter, type = "histogram")

sampled_invpdf_snps <- strsplit(rownames(taggingSNPs_v4_summary_invpdf), ":")

# Create a dataframe by applying a function to split_snps
sampled_snp_info <- do.call(rbind, lapply(sampled_invpdf_snps, function(snp) {
  chromosome <- snp[1]                 # Chromosome code (e.g., "I")
  start_range <- as.numeric(snp[2])     # Start of range (numeric position)
  
  # Apply the chromosome mapping
  #mapped_chromosome <- chromosome_map[[chromosome]]
  # End of range (start + 1)
  set_id <- paste0(chromosome, ":", start_range)  # Set ID
  
  # Return as a row in the dataframe
  return(data.frame(chromosome, set_id))
}))

write_tsv(sampled_snp_info, file = "/work/jb621/simHaplo/LD_mapping/ars_SS_u7o3_sims/StepWiseSuSiE/invpdf_sampling_causal/sampled_snp_info.tsv")
sampled_invpdf_snps <- as.data.frame(read_tsv(file = "/work/jb621/simHaplo/LD_mapping/ars_SS_u7o3_sims/StepWiseSuSiE/invpdf_sampling_causal/sampled_snp_info.tsv"))
