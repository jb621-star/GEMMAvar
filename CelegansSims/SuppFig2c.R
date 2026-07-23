# ===================================================
# Trying to map small h2 onto a more realistic scale
# ===================================================
# After your simulation loop, compute per-SNP h2
# which is what actually drives single-SNP power
r2_df <- data.frame(h2 = double(),
                    r2 = double())

per_snp_h2_list <- numeric(1000)

phen_means_list <- vector("list", n_sims)

for (s in 1:n_sims) {
  phen_means_list[[s]] <- read.table(
    paste0("/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/preWhiten/",
           "PhenotypeSimulated_singleSNP_h2005_UHetero1sc_sampledSNP_", s, ".txt")
  )$V1
}

for (s in 1:1000) {
  X  <- as.matrix(NORMAL.genomat.samp_mat[, s, drop = FALSE]) / 2
  Y  <- phen_means_list[[s]]   # already z-scaled
  
  # Regress Y on SNP, extract R-squared = per-SNP h2
  fit <- lm(Y ~ X)
  per_snp_h2_list[s] <- summary(fit)$r.squared
}

summary(per_snp_h2_list)
r2_df <- rbind(r2_df, data.frame(
  h2    = 0.05,
  r2    = per_snp_h2_list
))

save(r2_df, file = "/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/preWhiten/r2_vs_df.Rda")

# Plot
ggplot(r2_df, aes(x = r2, fill = h2, color = h2)) +
  geom_histogram(alpha = 0.3, position = "identity", bins = 50) +
  geom_vline(data = r2_df %>% group_by(h2) %>% summarise(mean = mean(r2)),
             aes(xintercept = mean, color = h2),
             linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  labs(
    title    = "Per-SNP R² distribution across heritability levels",
    subtitle = "Dashed lines show mean R² per h²",
    x        = "Per-SNP R² (causal variant)",
    y        = "Density",
    fill     = "Nominal h²",
    color    = "Nominal h²"
  ) +
  theme_bw() +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
        axis.title.x = element_text(size = 18),             # Increase x-axis label size
        axis.title.y = element_text(size = 18))             # Increase y-axis label size
