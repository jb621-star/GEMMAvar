library(tidyverse)
library(data.table)

load("/hpc/group/baughlab/jb621/GEMMA_2_testing/mouseTest/output/def_results_011326.Rda")
og_GWA_res <- read_tsv("mouse100_CD8MCH_lmm_demoP1Full.assoc.txt")

og_GWA_res <- og_GWA_res %>% 
  filter(rs %in% def_results$SNP)
plot(-log10(og_GWA_res$p_wald), -log10(results$p_wald))

head(og_GWA_res)
head(def_results)

combined <- og_GWA_res %>%
  left_join(def_results, by = c("rs" = "SNP"), suffix = c("_og", "_def")) %>% 
  mutate(log10_pval_Cpp = -log10(p_wald_og),
         log10_pval_R = -log10(p_wald_def))

CppvR <- ggplot(data = combined, aes(x = log10_pval_Cpp, y = log10_pval_R)) + geom_point(size = 2) +
  labs(y = "C++ Version of GEMMA -log10p", x = "R Versions of GEMMA -log10p", title = "Univariate Mouse CD8 p-value Comparison between GEMMA in C++ vs R") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated Negative Binomial with random strain effects, dispersion=0.04-0.05, logFC= 0.3 mean-var noise traits with 1000 different SNPs") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated Gamma(B = 0.25, s = 4) traits, heterogenous variance p-values using Garter Snake genotypes") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated h2 = 0.05 traits, imposed mean-var relationship p-values using Garter Snake genotypes") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated h2 = 0.025 traits, fitted mean-var relationship p-values using Garter Snake genotypes") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated h2 = 0.015 traits, fitted mean-var relationship p-values using Garter Snake genotypes") +
  geom_abline(slope = 1, color = "red") +  
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18))             # Increase y-axis label size
CppvR

load("UHetero_null_homog_results_060126.Rda")

NullHomo <- ggplot(data = power_results, aes(x = log10p_default, y = log10p_modified)) + geom_point(size = 2) +
  labs(y = "GEMMA -log10p", x = "GEMMAvar -log10p", title = "Simulated Null Normally Distributed traits with homogenous variances, tested on 1,000 different SNPs") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated Negative Binomial with random strain effects, dispersion=0.04-0.05, logFC= 0.3 mean-var noise traits with 1000 different SNPs") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated Gamma(B = 0.25, s = 4) traits, heterogenous variance p-values using Garter Snake genotypes") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated h2 = 0.05 traits, imposed mean-var relationship p-values using Garter Snake genotypes") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated h2 = 0.025 traits, fitted mean-var relationship p-values using Garter Snake genotypes") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated h2 = 0.015 traits, fitted mean-var relationship p-values using Garter Snake genotypes") +
  geom_abline(slope = 1, color = "red") +  
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18))             # Increase y-axis label size
NullHomo
