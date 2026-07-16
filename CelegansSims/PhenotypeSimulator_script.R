library(PhenotypeSimulator)
library(tidyverse)
library(vcfR)
library(ggplot2)
#library(genetics)
library(data.table)

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

  genetic_effects <- genBg_effect + genFixed$independent
  
  # ===== STEP 1: Generate base residual variances for non-outliers =====
  # =========>resid_var_base <- runif(n, 0, 1.5)
  
  # ===== STEP 2: Identify outlier individuals =====
  # Randomly select outlier indices (consistent across phenotype and variance)
  # ==========>outlier_idx <- sample(1:n, size = n_outliers, replace = FALSE)
  
  # ===== STEP 3: Modify variances for outliers =====
  # =========>resid_var_pred <- resid_var_base
  # =========>resid_var_pred[outlier_idx] <- resid_var_base[outlier_idx] * outlier_var_multiplier
  
  # 1. Heavy-tailed individual variances
  resid_var_base <- rlnorm(n, meanlog = -3, sdlog = 4)
  
  # Pareto variances
  #library(VGAM)
  #resid_var_base <- rpareto(n, shape = 0.8, scale = 0.5)  # shape < 1 = infinite variance
  
  # 2. Clip extreme values to stay numerically stable
  resid_var_base <- pmin(resid_var_base, quantile(resid_var_base, 0.98))
  resid_var_pred <- resid_var_base
  
  # 3. Draw noise preserving individual spread
  custom_noise <- rnorm(n, mean = 0, sd = sqrt(resid_var_pred))
  
  # Draw noise with heavy tails
  # custom_noise <- rt(n, df = 3) * sqrt(resid_var_base)
  
  noiseBg_scaled <- custom_noise * sqrt(noiseVar / var(custom_noise))
  
  Y_raw <- genetic_effects + noiseBg_scaled
  # Y_raw <- noiseBg_scaled + genBg_effect
  
  # diagnostics
  #plot(Y_raw, custom_noise, pch=20, xlab = "Simulated Phenotype", ylab = "Predicted Variance", main="Simulated Phenotypes vs Predicted variance")
  #plot(Y_raw, resid_var_pred, pch=20, xlab = "Real Slope Phenotype", ylab = "Predicted Variance", main="Real Phenotypes + Noise vs Predicted variance")
  
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
