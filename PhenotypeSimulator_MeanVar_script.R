library(PhenotypeSimulator)
library(tidyverse)
library(vcfR)
library(ggplot2)
#library(genetics)
library(data.table)

# Step 1: Load empirical data and fit mean-variance relationship
trait_file <- read_tsv("/hpc/group/baughlab/jb621/GEMMA_2_testing/slope100/full_traitData_slope.tsv", col_names = TRUE)
se2_file <- read_tsv("/hpc/group/baughlab/jb621/GEMMA_2_testing/slope100/Strain_TotalVar.tsv", col_names = FALSE)
trait_file$se2 <- se2_file$X1

hist(se2_file$X1)

# Fit log-variance ~ standardized mean relationship
mean_mean <- mean(trait_file$slope)
mean_sd   <- sd(trait_file$slope)
real_means_z <- (trait_file$slope - mean_mean) / mean_sd
eps_small <- 1e-8
mv_fit <- lm(log(trait_file$se2 + eps_small) ~ real_means_z)

plot(trait_file$slope, trait_file$se2, pch=20, main = "Mean-Residual Variance relationship in MIP Pilot data", xlab = "Average Slope values", ylab = "Total Residual Variance")
ord <- order(trait_file$slope)
lines(trait_file$slope[ord],
      exp(predict(mv_fit, newdata = data.frame(real_means_z = (trait_file$slope[ord] - mean_mean)/mean_sd))),
      col="red", lwd=2)
summary(mv_fit)   # check slope is negative as you observed

# Compare to the Z-scaled version of the data
# Step 1: Load empirical data and fit mean-variance relationship
# Load required data (your existing code)
load(file = "/work/jb621/simHaplo/LD_mapping/taggingSNPs_v4_geno_df.Rda")
#GarterfullGRM <- read.table("/hpc/group/baughlab/jb621/BulkGWA/MIPseqEvo/MIPseqEvo/20ds/gemma/grm_inbred/gemmaGRM.full.sXX.txt")
GarterfullGRM <- read.table("/work/jb621/simHaplo/LD_mapping/ars_SS_u7o3_sims_v2/ars_SS_u7o3_sims/StepWiseSuSiE/invpdf_sampling_causal/gemma/grm_inbred/gemmaGRM.full.sXX.txt", sep="\t", header=FALSE)
GarterfullGRM <- as.matrix(GarterfullGRM)
load(file = "/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/NORMgenomat.samp.Rda")

trait_file <- read_tsv("/hpc/group/baughlab/jb621/GEMMA_2_testing/slope100/full_traitData_slope.tsv", col_names = TRUE)
trait_file_zscale <- trait_file %>% 
  mutate(slope = as.numeric(scale(slope)))
se2_file <- read_tsv("/hpc/group/baughlab/jb621/GEMMA_2_testing/slope100/Strain_TotalVar_Zscale.tsv", col_names = TRUE)
trait_file_zscale$Variance <- se2_file$Variance

fit <- lm(log(Variance) ~ slope, data = trait_file_zscale)
confint(fit, level = 0.95) # Default level is 0.95, but can be specified

hist(trait_file_zscale$Variance)

# order the data
ord <- order(trait_file_zscale$slope)

summary(fit)

# predictions on the ordered data
pred <- predict(
  fit,
  newdata = data.frame(slope = trait_file_zscale$slope[ord]),
  interval = "confidence",
  level = 0.95
)

# if you fit log(Variance) ~ slope, back-transform with exp()
pred <- exp(pred)

# plot points
plot(trait_file_zscale$slope, trait_file_zscale$Variance, pch = 20,
     main = "Mean-Residual Variance relationship in Z-scaled MIP Pilot data",
     xlab = "Average Slope values", ylab = "Total Residual Variance")

# fitted line
lines(trait_file_zscale$slope[ord], pred[, "fit"], col = "red", lwd = 2)

# confidence interval (dashed)
lines(trait_file_zscale$slope[ord], pred[, "lwr"], col = "red", lwd = 1, lty = 2)
lines(trait_file_zscale$slope[ord], pred[, "upr"], col = "red", lwd = 1, lty = 2)

# Now apply the model to the Garter Snake standardized slope
Garter_slope <- read_tsv("/hpc/group/baughlab/jb621/GEMMA_2_testing/slope192/GarterMIP_slope_avg.tsv")
Garter_slope <- Garter_slope %>% 
  mutate(slope = as.numeric(scale(mean)))

ord <- order(Garter_slope$slope)

# predictions on the ordered data
pred <- predict(
  fit,
  newdata = data.frame(slope = Garter_slope$slope),
  interval = "confidence",
  level = 0.95
)

# if you fit log(Variance) ~ slope, back-transform with exp()
pred <- exp(pred)
Garter_predVar <- data.frame(pred)

# plot points
plot(Garter_slope$slope, Garter_predVar$fit, pch = 20,
     main = "Predicted Mean-Residual Variance relationship in Z-scaled Garter Snake data",
     xlab = "Average Slope values", ylab = "Predicted Total Residual Variance")

phen_means_raw_list <- vector("list", 1000)
phen_means_list     <- vector("list", 1000)
phen_vars_list      <- vector("list", 1000)
phen_noise_list     <- vector("list", 1000)
genetic_effects_list     <- vector("list", 1000)

genVar <- 0.05

# Null
# genVar <- 0.000001

# Non-null
#genVar <- 0.05
noiseVar <- 1 - genVar
h2s <- 0.99
phi <- 1
independent <- 1
n <- 192

# Compute eigendecomposition once outside the loop
eig <- eigen(GarterfullGRM, symmetric = TRUE)

# Pass to geneticBgEffects via a custom wrapper
# PhenotypeSimulator internally uses chol(), so bypass it by 
# simulating the background genetic effects manually

sim_bg_effects <- function(eig, n, genVar) {
  # Draw from MVN with covariance = kinship using eigendecomposition
  # equivalent to geneticBgEffects but without the chol() call
  sqrt_vals <- sqrt(pmax(eig$values, 0))
  z <- rnorm(n)
  effect <- eig$vectors %*% (sqrt_vals * z)
  # Scale to target variance
  effect <- as.vector(effect) * sqrt(genVar / var(effect))
  return(as.numeric(effect))
}

set.seed(43)

for (s in 1:1000) {
  #s <- 1
  # =========== GENETIC EFFECT ===============
  X_Garter_causal <- as.matrix(NORMAL.genomat.samp_mat[, s, drop = FALSE])/2
  genotype_vector <- X_Garter_causal
  
  #genFixed <- X_Garter_causal * 0.5
  genFixed <- geneticFixedEffects(
    X_causal = X_Garter_causal,
    P = 1, N = n, phenoID = "Trait",
    id_samples = rownames(X_Garter_causal),
    pTraitsAffected = 1, pIndependentGenetic = 1,
    pTraitIndependentGenetic = 1, keepSameIndependent = TRUE,
    distBeta = "norm", mBeta = 0.5, sdBeta = 0.1, verbose = FALSE
  )
  
  genBg_effect <- sim_bg_effects(eig, n = n, genVar = (1 - h2s) * genVar)
  
  genetic_effects <- genFixed$independent + genBg_effect
  
  # --- inside the loop, after genetic_effects computed ---
  
  # 2) rank-map genetic effects onto the empirical (z-scaled) slope distribution
  ranks <- rank(genetic_effects, ties.method = "first")
  probs <- (ranks - 0.5) / n
  mapped_slopes <- as.numeric(quantile(trait_file_zscale$slope, probs = probs, type = 8))
  
  # 3) predict variance from the fitted mean-variance model at these mapped positions
  pred <- predict(fit, newdata = data.frame(slope = mapped_slopes),
                  interval = "confidence", level = 0.95)
  resid_var_pred <- exp(pred[, "fit"])   # back-transform since fit was on log(Variance)
  
  # 4) guard against pathological values
  if (any(is.na(resid_var_pred) | resid_var_pred <= 0)) {
    fallback <- min(resid_var_pred[resid_var_pred > 0], na.rm = TRUE)
    resid_var_pred[is.na(resid_var_pred) | resid_var_pred <= 0] <- fallback
  }
  
  
  # ===== STEP 1: Generate base residual variances for non-outliers =====
  # =========>resid_var_base <- runif(n, 0, 1.5)
  
  # ===== STEP 2: Identify outlier individuals =====
  # Randomly select outlier indices (consistent across phenotype and variance)
  # ==========>outlier_idx <- sample(1:n, size = n_outliers, replace = FALSE)
  
  # ===== STEP 3: Modify variances for outliers =====
  # =========>resid_var_pred <- resid_var_base
  # =========>resid_var_pred[outlier_idx] <- resid_var_base[outlier_idx] * outlier_var_multiplier
  
  # noiseBg_scaled <- custom_noise * (noiseVar / mean(custom_noise))
  
  # 5) draw heteroskedastic noise (use these variances only)
  # ==============> custom_noise <- rnorm(n, mean = 0, sd = sqrt(resid_var_pred))
  #custom_noise <- rnorm(n, mean = 0, sd = runif(1, min = 0, max = 1.5))
  # scale using realized variance, not theoretical
  # ==============> noiseBg_scaled <- custom_noise * sqrt(noiseVar / var(custom_noise))
  
  # # noiseBg <- noiseBgEffects(n, P = 1, mean = 0, sd = sqrt(Garter_predVar$fit), shared = FALSE, independent = TRUE)
  # # noiseBg_scaled <- rescaleVariance(noiseBg$independent, independent * phi * noiseVar)
  # 
  # # draw heteroskedastic noise directly
  # # predicted from fit using some proxy mean (e.g., mapped_means or Garter_predVar)
  # # predict residual variances using mapped_means
  # # 
  
  # Heavy-tailed individual variances
  # resid_var_base <- rlnorm(n, meanlog = -3, sdlog = 4)
  
  
  # Draw noise preserving individual spread
  custom_noise <- rnorm(n, mean = 0, sd = sqrt(resid_var_pred))
  
  noiseBg_scaled <- custom_noise * sqrt(noiseVar / var(custom_noise))
  
  # custom_noise_centered <- custom_noise - mean(custom_noise)
  # noiseBg_scaled <- custom_noise_centered * sqrt(noiseVar / var(custom_noise_centered))
  
  # # 
  # # # rescale to hit desired global noise variance
  # # current_noise_var <- var(custom_noiseBg)
  # # scale_factor <- sqrt(noiseVar / current_noise_var)
  # # noiseBg_scaled <- custom_noiseBg * scale_factor
  # # Y_raw <- genBg$independent[,1] + 
  # #   genFixed$independent[,1] +
  # #   custom_noise
  # 
  # # Y_raw <- genBg_scaled$component[, 1] +
  # #      genFixed_scaled$component[, 1] +
  # #   custom_noise
  # # Combine to produce total phenotype
  # #Y_slope <- Garter_slope$slope + custom_noiseBg
  
  Y_raw <- genetic_effects + noiseBg_scaled
  #Y_raw <- noiseBg_scaled + genBg_effect
  
  # diagnostics
  plot(Y_raw, custom_noise, pch=20, xlab = "Simulated Value", ylab = "Noise Term (e ~ N(0, Predicted Variance))", main="Simulated Phenotypes vs Noise Term")
  plot(Y_raw, resid_var_pred, pch=20, xlab = "Simulated Trait Value", ylab = "Predicted Variance", main="Simulated Phenotype vs Predicted variance")
  
  # Store results
  phen_means_raw_list[[s]] <- Y_raw
  phen_means_list[[s]] <- scale(Y_raw)[,1]
  phen_vars_list[[s]]      <- resid_var_pred
  phen_noise_list[[s]]     <- custom_noise
}

hist(phen_vars_list[[1]])

# Diagnostics
cat("group counts:", table(geno_vec), "\n")
cat("mean resid_var by group:", tapply(resid_var, geno_vec, mean), "\n")
cat("observed noise var by group:", tapply(custom_noiseBg, geno_vec, var), "\n")
cat("genetic var:", var(gvec), "noise var:", var(custom_noiseBg), "SNR:", var(gvec)/var(custom_noiseBg), "\n")

print(summary(lm(Y_raw ~ geno_vec)))
print(summary(lm(Y_raw ~ geno_vec, weights = 1/resid_var)))
# weighted version fits better...so what happened???

# Diagnostic plots
plot(phen_means_list[[1]], phen_vars_list[[1]],  xlab = "Mean phenotypes", ylab = "Residual Variance", main = "PhenotypeSimulator with null random residual variance log-normal(m = 0 ,var = 1.5), Iteration 1")
var(phen_vars_list[[1]])
hist(phen_means_list[[1]])

par(mfrow = c(1, 3))
hist(se2_file$X1, main = "Input Residual Variances", xlab = "Variance")
hist(phen_means_list[[1]], main = "Average Trait Values (Sim 1)", xlab = "Trait Value")
hist(phen_vars_list[[1]], main = "Empirical Within-Strain Variances (Sim 1)", xlab = "Variance")
par(mfrow = c(1, 1))
boxplot(phen_means_list[[1]] ~ NORMAL.genomat.samp_mat[,1], xlab = "Allele at Causal SNP 1", main = "Phenotype by Genotype (Sim 1)", ylab = "Scaled Phenotype")
#boxplot(mapped_genetic_effects ~ NORMAL.genomat.samp_mat[,1], xlab = "Allele at Causal SNP 1", main = "Phenotype by Genotype (Sim 1)", ylab = "Scaled Phenotype")

# Write Null Phentoype files
for (s in 1:1000) {
  write.table(phen_means_list[[s]], 
              file = paste0("/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/preWhiten/PhenotypeSimulated_singleSNP_null_rlnormHETERO_mneg3sd4_0602_sampledSNP_", s, ".txt"), 
              col.names = FALSE, row.names = FALSE)
  write.table(phen_vars_list[[s]],
              file = paste0("/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/preWhiten/PhenotypeSimulated_singleSNP_residSE2_null_rlnormHETERO_mneg3sd4_0602_sampledSNP_", s, ".txt"), 
              col.names = FALSE, row.names = FALSE)
}
