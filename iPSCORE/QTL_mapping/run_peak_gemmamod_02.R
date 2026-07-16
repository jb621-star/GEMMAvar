## ============================================================
## 02_run_peak_gemmamod_caQTL.R
## Called by SLURM array — runs GEMMAmod for all SNPs of one ATAC-seq peak
## Usage: Rscript 02_run_peak_gemmamod_caQTL.R <peak_id>
## ============================================================

.libPaths(c("/home/jb621/R/x86_64-pc-linux-gnu-library/4.2",
            "/usr/local/lib/R/site-library",
            "/usr/local/lib/R/library",
            .libPaths()))

library(limma)
library(edgeR)
library(tidyverse)
library(data.table)
library(parallel)

# ============================================================================
# Core REML Likelihood Functions
# ============================================================================

likelihood_derivative1_lambda <- function(lambda, eigenVals, Y, W) {
  n <- nrow(Y)
  c <- ncol(W)
  D <- lambda * eigenVals + 1.0
  W_scaled <- W / sqrt(D)
  WtW_inv <- solve(crossprod(W_scaled))
  P_diag <- 1/D - rowSums((W/D) * (W %*% WtW_inv %*% t(W/D)))
  trace_term <- sum(eigenVals * P_diag)
  quad_term <- sum((eigenVals * P_diag) * (Y^2))
  deriv <- -0.5 * trace_term + 0.5 * quad_term / sum(P_diag * Y^2) * (n - c)
  return(deriv)
}

likelihood_derivative2_lambda <- function(lambda, eigenVals, Y, W) {
  h <- 1e-5
  deriv_plus  <- likelihood_derivative1_lambda(lambda + h, eigenVals, Y, W)
  deriv_minus <- likelihood_derivative1_lambda(lambda - h, eigenVals, Y, W)
  return((deriv_plus - deriv_minus) / (2 * h))
}

likelihood_lambda <- function(lambda, eigenVals, Y, W) {
  n <- nrow(Y)
  c <- ncol(W)
  D <- lambda * eigenVals + 1.0
  logdet_D <- sum(log(D))
  W_scaled <- W / sqrt(D)
  WtDW <- crossprod(W_scaled)
  logdet_WtDW <- determinant(WtDW, logarithm = TRUE)$modulus[1]
  Y_D <- Y / D
  WtDW_inv <- solve(WtDW)
  Y_P_Y <- sum(Y * Y_D) - t(crossprod(W, Y_D)) %*% WtDW_inv %*% crossprod(W, Y_D)
  loglik <- -0.5 * (logdet_D + logdet_WtDW + (n - c) * log(Y_P_Y[1]))
  return(as.numeric(loglik))
}

# ============================================================================
# Lambda Optimization
# ============================================================================

calc_lambda_restricted <- function(eigenVals, Y, W, grid = FALSE) {
  if (grid) {
    step <- 1.0
    lambda_pow_low  <- -5.0
    lambda_pow_high <-  5.0
    lambda_grid <- seq(lambda_pow_low, lambda_pow_high - step, by = step)
    roots <- c(10^lambda_pow_low, 10^lambda_pow_high)
    likelihood_prev <- likelihood_derivative1_lambda(10^lambda_pow_low, eigenVals, Y, W)
    for (i in seq_along(lambda_grid)) {
      lambda0 <- 10^lambda_grid[i]
      lambda1 <- 10^(lambda_grid[i] + step)
      if (i > 1) {
        likelihood_lambda0 <- likelihood_lambda1
      } else {
        likelihood_lambda0 <- likelihood_prev
      }
      likelihood_lambda1 <- likelihood_derivative1_lambda(lambda1, eigenVals, Y, W)
      if (sign(likelihood_lambda0) * sign(likelihood_lambda1) < 0) {
        lambda_min <- tryCatch({
          uniroot(f = function(l) likelihood_derivative1_lambda(l, eigenVals, Y, W),
                  lower = lambda0, upper = lambda1, tol = 0.1)$root
        }, error = function(e) { (lambda0 + lambda1) / 2 })
        lambda_min <- tryCatch({
          for (iter in 1:10) {
            f_val   <- likelihood_derivative1_lambda(lambda_min, eigenVals, Y, W)
            f_prime <- likelihood_derivative2_lambda(lambda_min, eigenVals, Y, W)
            if (abs(f_prime) < 1e-10) break
            lambda_new <- lambda_min - f_val / f_prime
            if (abs(lambda_new - lambda_min) < 1e-5) break
            lambda_min <- lambda_new
          }
          lambda_min
        }, error = function(e) { lambda_min })
        roots <- c(roots, lambda_min)
      }
    }
    likelihoods <- sapply(roots, function(lam) likelihood_lambda(lam, eigenVals, Y, W))
    return(roots[which.max(likelihoods)])
  } else {
    opt_result <- tryCatch({
      optimize(f = function(l) -likelihood_lambda(l, eigenVals, Y, W),
               lower = 1e-5, upper = 1e5, tol = 1e-4)
    }, error = function(e) { list(minimum = 1.0) })
    return(opt_result$minimum)
  }
}

# ============================================================================
# Parameter Estimation
# ============================================================================

calc_beta_vg_ve_restricted <- function(eigenVals, W, X, lambda, Y) {
  n <- nrow(Y)
  c <- ncol(W)
  D <- lambda * eigenVals + 1.0
  WX <- cbind(W, X)
  WX_D <- WX / sqrt(D)
  WXtDWX <- crossprod(WX_D)
  WXtDWX_inv <- solve(WXtDWX)
  WXtDY <- crossprod(WX / D, Y)
  beta_vec <- WXtDWX_inv %*% WXtDY
  beta <- beta_vec[nrow(beta_vec), 1]
  Y_D <- Y / D
  fitted <- WX %*% beta_vec
  fitted_D <- fitted / D
  Y_P_Y <- sum(Y * Y_D) - sum(fitted * fitted_D)
  tau <- (n - c - 1) / Y_P_Y
  se_beta <- sqrt(WXtDWX_inv[nrow(WXtDWX_inv), ncol(WXtDWX_inv)] / tau)
  return(list(beta    = as.numeric(beta),
              beta_vec = beta_vec,
              se_beta = as.numeric(se_beta),
              tau     = as.numeric(tau)))
}

# ============================================================================
# Main pyGEMMA Function
# ============================================================================

pygemma <- function(Y, X, W, K, Z = NULL, snps = NULL,
                    verbose = 0, disable_checks = TRUE,
                    grid = FALSE, eigen = TRUE, nproc = 1) {
  Y <- as.matrix(Y); X <- as.matrix(X)
  W <- as.matrix(W); K <- as.matrix(K)
  if (ncol(Y) != 1) Y <- matrix(Y, ncol = 1)
  if (!is.null(Z)) { Z <- as.matrix(Z); K <- Z %*% K %*% t(Z) }
  n <- nrow(Y); m <- ncol(X); c <- ncol(W)
  if (eigen) {
    eigen_decomp <- eigen(K, symmetric = TRUE)
    eigenVals <- pmax(0, eigen_decomp$values)
    U <- eigen_decomp$vectors
    X <- t(U) %*% X; Y <- t(U) %*% Y; W <- t(U) %*% W
  } else {
    eigenVals <- pmax(0, K)
  }
  if (!disable_checks) {
    if (any(is.nan(X)) || any(is.nan(Y)) || any(is.nan(W))) stop("NaN values present")
  }
  process_snp <- function(g) {
    x_g <- X[, g, drop = FALSE]
    tryCatch({
      lambda_restricted <- calc_lambda_restricted(eigenVals, Y, cbind(W, x_g), grid = grid)
      params <- calc_beta_vg_ve_restricted(eigenVals, W, x_g, lambda_restricted, Y)
      F_wald <- (params$beta / params$se_beta)^2
      p_wald <- pf(F_wald, df1 = 1, df2 = n - c - 1, lower.tail = FALSE)
      data.frame(beta = params$beta, se_beta = params$se_beta,
                 tau = params$tau, lambda = lambda_restricted,
                 F_wald = F_wald, p_wald = p_wald,
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      if (verbose > 1) cat(sprintf("Error SNP %d: %s\n", g, e$message))
      data.frame(beta = NA_real_, se_beta = NA_real_, tau = NA_real_,
                 lambda = NA_real_, F_wald = NA_real_, p_wald = NA_real_,
                 stringsAsFactors = FALSE)
    })
  }
  if (nproc > 1 && m > nproc) {
    cl <- makeCluster(min(nproc, m))
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(cl, c("eigenVals","Y","W","X","n","c","grid",
                        "calc_lambda_restricted","calc_beta_vg_ve_restricted",
                        "likelihood_derivative1_lambda","likelihood_derivative2_lambda",
                        "likelihood_lambda"), envir = environment())
    results_list <- parLapply(cl, 1:m, process_snp)
  } else {
    results_list <- lapply(1:m, process_snp)
  }
  results_df <- do.call(rbind, results_list)
  if (!is.null(snps)) results_df$SNP <- snps
  return(results_df)
}

# ============================================================================
# Main script
# ============================================================================

setwd("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/")

args    <- commandArgs(trailingOnly = TRUE)
#peak_id <- peak_list[1]
peak_id <- args[1]

if (is.na(peak_id) || peak_id == "") {
  stop("No peak_id supplied. Usage: Rscript 02_run_peak_gemmamod_caQTL.R <peak_id>")
}

input_file  <- paste0("caqtl_gene_inputs/",        peak_id, ".rds")
output_file <- paste0("caqtl_gemmamod_outputs/",   peak_id, "NOSHRINKAGE.rds")

dir.create("caqtl_gemmamod_outputs", showWarnings = FALSE)

if (file.exists(output_file)) {
  cat("Output already exists, skipping:", peak_id, "\n")
  quit(status = 0)
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Load shared inputs
W        <- readRDS("caqtl_gene_inputs/W.rds")
K        <- readRDS("caqtl_gene_inputs/K.rds")
dge_atac <- readRDS("caqtl_gene_inputs/dge_atac.rds")  # DGEList for ATAC peaks

#dge_atac_samples <- colnames(dge_atac)

# dge_samples <- colnames(dge_atac)
# 
# cat("Samples in ATAC not in genotypes:", 
#     sum(!dge_samples %in% geno_samples), "\n")
# cat("Samples in genotypes not in ATAC:", 
#     sum(!geno_samples %in% dge_samples), "\n")
# cat("Samples in common:", 
#     sum(dge_samples %in% geno_samples), "\n")

# Load peak-specific inputs
inp     <- readRDS(input_file)
geno    <- inp$geno
snp_ids <- inp$snp_ids
n_snps  <- length(snp_ids)

#geno_samples <- rownames(geno)
#K_samples <- colnames(K)

# atac_to_geno <- meta_map$ATAC_UUID[match(geno_samples, meta_map$WGS_UUID)]
# dge_atac <- dge_atac[, atac_to_geno]
# colnames(dge_atac)
# saveRDS(dge_atac, "caqtl_gene_inputs/dge_atac.rds")

cat(sprintf("Peak: %s | SNPs: %d\n", peak_id, n_snps))

# Pre-allocate results
mod_results <- data.frame(
  SNP_ID     = character(n_snps),
  Element_ID = character(n_snps),
  beta       = numeric(n_snps),
  se_beta    = numeric(n_snps),
  tau        = numeric(n_snps),
  lambda     = numeric(n_snps),
  F_wald     = numeric(n_snps),
  p_wald     = numeric(n_snps),
  stringsAsFactors = FALSE
)

# ---- Per-SNP loop ----
for (i in seq_len(n_snps)) {
  
  snp_id      <- snp_ids[i]
  causal_geno <- geno[, snp_id]
  
  if (any(is.na(causal_geno))) {
    causal_geno[is.na(causal_geno)] <- mean(causal_geno, na.rm = TRUE)
  }
  
  design_voom <- model.matrix(~ causal_geno)
  y           <- voom(dge_atac, design = design_voom, plot = FALSE)
  
  expr_peak    <- y$E[peak_id, , drop = FALSE]                                     # 1 x n
  weights_peak <- y$weights[match(peak_id, rownames(dge_atac)), , drop = FALSE]    # 1 x n
  
  # Shrink log-variances
  se_pheno      <- as.matrix(1 / weights_peak)
  log_s2        <- log(se_pheno)
  mu            <- mean(log_s2)
  lambda_shrink <- 0
  log_s2_shr    <- (1 - lambda_shrink) * log_s2 + lambda_shrink * mu
  s2_shr        <- exp(log_s2_shr)
  
  w   <- 1 / sqrt(s2_shr)
  Ww  <- diag(as.vector(w))
  
  K_whitened     <- Ww %*% as.matrix(K) %*% Ww
  Y_w            <- Ww %*% t(expr_peak)
  X_w            <- Ww %*% as.matrix(causal_geno, ncol = 1)
  W_fe_w         <- Ww %*% W
  K_whitened_reg <- K_whitened + diag(1e-4, nrow(K_whitened))
  
  result <- tryCatch(
    pygemma(Y = Y_w, X = X_w, W = W_fe_w, K = K_whitened_reg,
            snps = snp_id, verbose = 0, grid = FALSE, nproc = 1),
    error = function(e) {
      cat(sprintf("  ERROR SNP %s: %s\n", snp_id, conditionMessage(e)))
      NULL
    }
  )
  
  mod_results$SNP_ID[i]     <- snp_id
  mod_results$Element_ID[i] <- peak_id
  
  if (!is.null(result)) {
    mod_results$beta[i]    <- result$beta
    mod_results$se_beta[i] <- result$se_beta
    mod_results$tau[i]     <- result$tau
    mod_results$lambda[i]  <- result$lambda
    mod_results$F_wald[i]  <- result$F_wald
    mod_results$p_wald[i]  <- result$p_wald
  } else {
    mod_results$beta[i]    <- NA
    mod_results$se_beta[i] <- NA
    mod_results$tau[i]     <- NA
    mod_results$lambda[i]  <- NA
    mod_results$F_wald[i]  <- NA
    mod_results$p_wald[i]  <- NA
  }
}

saveRDS(mod_results, output_file)
cat(sprintf("Done: %s — %d SNPs tested, %d errors\n",
            peak_id, n_snps, sum(is.na(mod_results$beta))))