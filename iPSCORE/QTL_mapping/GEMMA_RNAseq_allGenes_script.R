.libPaths(c("/home/jb621/R/x86_64-pc-linux-gnu-library/4.2","/usr/local/lib/R/site-library","/usr/local/lib/R/library",.libPaths()) ) 
library(limma)
library(edgeR)
library(tidyverse)
#library(Matrix)
library(data.table)
library(parallel)

setwd("~/QTL_mapping/")

# ============================================================================
# Core REML Likelihood Functions
# ============================================================================

#' First derivative of REML likelihood with respect to lambda
#'
#' @param lambda Variance component ratio
#' @param eigenVals Eigenvalues of kinship matrix
#' @param Y Transformed phenotype vector
#' @param W Transformed covariate matrix
#' @return First derivative value
likelihood_derivative1_lambda <- function(lambda, eigenVals, Y, W) {
  n <- nrow(Y)
  c <- ncol(W)
  
  # Compute D = lambda * eigenVals + 1
  D <- lambda * eigenVals + 1.0
  
  # Compute P = I - W(W'D^-1W)^-1W'D^-1
  W_scaled <- W / sqrt(D)
  WtW_inv <- solve(crossprod(W_scaled))
  
  # P_diag represents diagonal of P
  P_diag <- 1/D - rowSums((W/D) * (W %*% WtW_inv %*% t(W/D)))
  
  # Compute derivative
  trace_term <- sum(eigenVals * P_diag)
  quad_term <- sum((eigenVals * P_diag) * (Y^2))
  
  deriv <- -0.5 * trace_term + 0.5 * quad_term / sum(P_diag * Y^2) * (n - c)
  
  return(deriv)
}

#' Second derivative of REML likelihood with respect to lambda
#'
#' @param lambda Variance component ratio
#' @param eigenVals Eigenvalues of kinship matrix
#' @param Y Transformed phenotype vector
#' @param W Transformed covariate matrix
#' @return Second derivative value
likelihood_derivative2_lambda <- function(lambda, eigenVals, Y, W) {
  n <- nrow(Y)
  c <- ncol(W)
  
  D <- lambda * eigenVals + 1.0
  
  # More complex computation for second derivative
  # This is a simplified version - full implementation would need careful derivation
  W_scaled <- W / sqrt(D)
  WtW_inv <- solve(crossprod(W_scaled))
  P_diag <- 1/D - rowSums((W/D) * (W %*% WtW_inv %*% t(W/D)))
  
  # Approximate second derivative
  h <- 1e-5
  deriv_plus <- likelihood_derivative1_lambda(lambda + h, eigenVals, Y, W)
  deriv_minus <- likelihood_derivative1_lambda(lambda - h, eigenVals, Y, W)
  
  return((deriv_plus - deriv_minus) / (2 * h))
}

#' REML likelihood function
#'
#' @param lambda Variance component ratio
#' @param eigenVals Eigenvalues of kinship matrix
#' @param Y Transformed phenotype vector
#' @param W Transformed covariate matrix
#' @return REML log-likelihood value
likelihood_lambda <- function(lambda, eigenVals, Y, W) {
  n <- nrow(Y)
  c <- ncol(W)
  
  D <- lambda * eigenVals + 1.0
  
  # Log determinant term
  logdet_D <- sum(log(D))
  
  # Compute W'D^-1W
  W_scaled <- W / sqrt(D)
  WtDW <- crossprod(W_scaled)
  logdet_WtDW <- determinant(WtDW, logarithm = TRUE)$modulus[1]
  
  # Compute Y'PY where P = D^-1 - D^-1W(W'D^-1W)^-1W'D^-1
  Y_D <- Y / D
  WtDW_inv <- solve(WtDW)
  Y_P_Y <- sum(Y * Y_D) - t(crossprod(W, Y_D)) %*% WtDW_inv %*% crossprod(W, Y_D)
  
  # REML log-likelihood
  loglik <- -0.5 * (logdet_D + logdet_WtDW + (n - c) * log(Y_P_Y[1]))
  
  return(as.numeric(loglik))
}

# ============================================================================
# Lambda Optimization
# ============================================================================

#' Calculate optimal lambda using grid search and optimization
#'
#' @param eigenVals Eigenvalues of kinship matrix
#' @param Y Transformed phenotype vector
#' @param W Transformed covariate matrix (including SNP for restricted model)
#' @param grid Logical, use grid search (slower but more robust)
#' @return Optimal lambda value
calc_lambda_restricted <- function(eigenVals, Y, W, grid = FALSE) {
  
  if (grid) {
    # Grid search approach (more robust)
    step <- 1.0
    lambda_pow_low <- -5.0
    lambda_pow_high <- 5.0
    
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
      
      # Check for sign change
      if (sign(likelihood_lambda0) * sign(likelihood_lambda1) < 0) {
        # Use Brent's method to find root
        lambda_min <- tryCatch({
          uniroot(
            f = function(l) likelihood_derivative1_lambda(l, eigenVals, Y, W),
            lower = lambda0,
            upper = lambda1,
            tol = 0.1
          )$root
        }, error = function(e) {
          (lambda0 + lambda1) / 2
        })
        
        # Refine with Newton's method
        lambda_min <- tryCatch({
          for (iter in 1:10) {
            f_val <- likelihood_derivative1_lambda(lambda_min, eigenVals, Y, W)
            f_prime <- likelihood_derivative2_lambda(lambda_min, eigenVals, Y, W)
            
            if (abs(f_prime) < 1e-10) break
            
            lambda_new <- lambda_min - f_val / f_prime
            
            if (abs(lambda_new - lambda_min) < 1e-5) break
            
            lambda_min <- lambda_new
          }
          lambda_min
        }, error = function(e) {
          lambda_min
        })
        
        roots <- c(roots, lambda_min)
      }
    }
    
    # Evaluate likelihood at all roots
    likelihoods <- sapply(roots, function(lam) {
      likelihood_lambda(lam, eigenVals, Y, W)
    })
    
    return(roots[which.max(likelihoods)])
    
  } else {
    # Fast optimization approach
    opt_result <- tryCatch({
      optimize(
        f = function(l) -likelihood_lambda(l, eigenVals, Y, W),
        lower = 1e-5,
        upper = 1e5,
        tol = 1e-4
      )
    }, error = function(e) {
      list(minimum = 1.0)
    })
    
    return(opt_result$minimum)
  }
}

# ============================================================================
# Parameter Estimation
# ============================================================================

#' Calculate beta, standard error, and variance components (REML)
#'
#' @param eigenVals Eigenvalues of kinship matrix
#' @param W Transformed covariate matrix
#' @param X Transformed SNP genotypes (single SNP)
#' @param lambda Variance component ratio
#' @param Y Transformed phenotype vector
#' @return List containing beta, beta_vec, se_beta, tau
calc_beta_vg_ve_restricted <- function(eigenVals, W, X, lambda, Y) {
  n <- nrow(Y)
  c <- ncol(W)
  
  # Compute D = lambda * eigenVals + 1
  D <- lambda * eigenVals + 1.0
  
  # Form augmented covariate matrix [W, X]
  WX <- cbind(W, X)
  
  # Compute (W|X)'D^-1(W|X)
  WX_D <- WX / sqrt(D)
  WXtDWX <- crossprod(WX_D)
  WXtDWX_inv <- solve(WXtDWX)
  
  # Compute beta vector: (W|X)'D^-1Y
  WXtDY <- crossprod(WX / D, Y)
  beta_vec <- WXtDWX_inv %*% WXtDY
  
  # Extract SNP effect (last element)
  beta <- beta_vec[nrow(beta_vec), 1]
  
  # Compute residual variance (REML)
  # Y'PY where P = D^-1 - D^-1(W|X)[(W|X)'D^-1(W|X)]^-1(W|X)'D^-1
  Y_D <- Y / D
  fitted <- WX %*% beta_vec
  fitted_D <- fitted / D
  Y_P_Y <- sum(Y * Y_D) - sum(fitted * fitted_D)
  
  tau <- (n - c - 1) / Y_P_Y
  
  # Standard error of beta (last diagonal element)
  se_beta <- sqrt(WXtDWX_inv[nrow(WXtDWX_inv), ncol(WXtDWX_inv)] / tau)
  
  return(list(
    beta = as.numeric(beta),
    beta_vec = beta_vec,
    se_beta = as.numeric(se_beta),
    tau = as.numeric(tau)
  ))
}

# ============================================================================
# Main pyGEMMA Function
# ============================================================================

#' Perform GEMMA analysis using REML
#'
#' @param Y Phenotype matrix (n x 1)
#' @param X Genotype matrix (n x m), where m is number of SNPs
#' @param W Covariate matrix (n x c)
#' @param K Genetic relatedness matrix (GRM) (n x n)
#' @param Z Optional random effect matrix (default NULL)
#' @param snps Vector of SNP names (default NULL)
#' @param verbose Verbosity level (0 = silent, 1 = progress, 2 = detailed)
#' @param disable_checks Disable NaN checks for speed (default TRUE)
#' @param grid Use grid search for lambda (slower but more robust, default FALSE)
#' @param eigen Use eigendecomposition (default TRUE)
#' @param nproc Number of parallel processes (default 1)
#' @return Data frame with results
pygemma <- function(Y, X, W, K, Z = NULL, snps = NULL, 
                    verbose = 0, disable_checks = TRUE, 
                    grid = FALSE, eigen = TRUE, nproc = 1) {
  
  # Convert to matrices and ensure proper types
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  W <- as.matrix(W)
  K <- as.matrix(K)
  
  if (ncol(Y) != 1) {
    Y <- matrix(Y, ncol = 1)
  }
  
  # Apply Z transformation if provided
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    K <- Z %*% K %*% t(Z)
  }
  
  n <- nrow(Y)
  m <- ncol(X)
  c <- ncol(W)
  
  if (verbose > 0) {
    cat("Running pyGEMMA with REML...\n")
    cat(sprintf("  Samples: %d, SNPs: %d, Covariates: %d\n", n, m, c))
  }
  
  # Eigendecomposition
  if (verbose > 0) cat("Performing eigendecomposition...\n")
  
  if (eigen) {
    eigen_decomp <- eigen(K, symmetric = TRUE)
    eigenVals <- pmax(0, eigen_decomp$values)  # Set negative eigenvalues to 0
    U <- eigen_decomp$vectors
    
    if (verbose > 0) cat("  Eigendecomposition complete\n")
    
    # Transform data
    if (verbose > 0) cat("Transforming data...\n")
    X <- t(U) %*% X
    Y <- t(U) %*% Y
    W <- t(U) %*% W
    
  } else {
    # Use K directly (set negative values to 0)
    K <- pmax(0, K)
    eigenVals <- K
  }
  
  # Check for NaN values
  if (!disable_checks) {
    if (any(is.nan(X)) || any(is.nan(Y)) || any(is.nan(W))) {
      stop("NaN values present in data")
    }
  }
  
  if (verbose > 0) cat(sprintf("Testing %d SNPs...\n", m))
  
  # Define function to process each SNP
  process_snp <- function(g) {
    x_g <- X[, g, drop = FALSE]
    
    tryCatch({
      # Calculate lambda for restricted model
      lambda_restricted <- calc_lambda_restricted(eigenVals, Y, cbind(W, x_g), grid = grid)
      
      # Calculate beta, se, tau
      params <- calc_beta_vg_ve_restricted(eigenVals, W, x_g, lambda_restricted, Y)
      
      # Wald test
      F_wald <- (params$beta / params$se_beta)^2
      p_wald <- pf(F_wald, df1 = 1, df2 = n - c - 1, lower.tail = FALSE)
      
      data.frame(
        beta = params$beta,
        se_beta = params$se_beta,
        tau = params$tau,
        lambda = lambda_restricted,
        F_wald = F_wald,
        p_wald = p_wald,
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      if (verbose > 1) cat(sprintf("Error processing SNP %d: %s\n", g, e$message))
      data.frame(
        beta = NA_real_,
        se_beta = NA_real_,
        tau = NA_real_,
        lambda = NA_real_,
        F_wald = NA_real_,
        p_wald = NA_real_,
        stringsAsFactors = FALSE
      )
    })
  }
  
  # Process SNPs (with or without parallelization)
  if (nproc > 1 && m > nproc) {
    if (verbose > 0) cat(sprintf("Using %d parallel processes...\n", nproc))
    
    # Setup cluster
    cl <- makeCluster(min(nproc, m))
    on.exit(stopCluster(cl), add = TRUE)
    
    # Export necessary objects and functions
    clusterExport(cl, c("eigenVals", "Y", "W", "X", "n", "c", "grid",
                        "calc_lambda_restricted", "calc_beta_vg_ve_restricted",
                        "likelihood_derivative1_lambda", "likelihood_derivative2_lambda",
                        "likelihood_lambda"),
                  envir = environment())
    
    # Run in parallel with progress
    if (verbose > 0) {
      results_list <- pbapply::pblapply(1:m, process_snp, cl = cl)
    } else {
      results_list <- parLapply(cl, 1:m, process_snp)
    }
    
  } else {
    # Sequential processing
    if (verbose > 0 && requireNamespace("pbapply", quietly = TRUE)) {
      results_list <- pbapply::pblapply(1:m, process_snp)
    } else {
      results_list <- lapply(1:m, process_snp)
    }
  }
  
  # Combine results
  results_df <- do.call(rbind, results_list)
  
  # Add SNP names if provided
  if (!is.null(snps)) {
    results_df$SNP <- snps
  }
  
  if (verbose > 0) cat("Analysis complete!\n")
  
  return(results_df)
}

# ============================================================================
# Helper function to load and prepare data
# ============================================================================

#' Prepare data for pyGEMMA analysis
#'
#' @param pheno Phenotype vector or data frame
#' @param geno Genotype matrix
#' @param covars Covariate matrix (intercept will be added if not present)
#' @param kinship Kinship matrix
#' @return List with prepared Y, X, W, K matrices
prepare_gemma_data <- function(pheno, geno, covars = NULL, kinship) {
  
  # Prepare phenotype
  Y <- as.matrix(pheno)
  if (ncol(Y) != 1) Y <- matrix(Y, ncol = 1)
  
  # Prepare genotypes
  X <- as.matrix(geno)
  
  # Prepare covariates (add intercept if needed)
  if (is.null(covars)) {
    W <- matrix(1, nrow = nrow(Y), ncol = 1)
  } else {
    W <- as.matrix(covars)
    # Check if intercept exists (column of all 1s)
    has_intercept <- any(apply(W, 2, function(col) all(col == col[1])))
    if (!has_intercept) {
      W <- cbind(1, W)
    }
  }
  
  # Prepare kinship
  K <- as.matrix(kinship)
  
  # Check dimensions
  n <- nrow(Y)
  if (nrow(X) != n || nrow(W) != n || nrow(K) != n || ncol(K) != n) {
    stop("Dimension mismatch in input matrices")
  }
  
  list(Y = Y, X = X, W = W, K = K)
}

# Phenotypes
# # Normalized + filtered + INT-transformed matrix (for eQTL input)
# edgeR_expr <- read.table(file = "ppc_rna_int_normalized.txt",
#             sep       = "\t",
#             header = TRUE,
#             row.names = 1,
#             check.names = FALSE)   # col.names = NA writes row names correctly

# VOOM CPM (post-filter) 
# it might be incorrect to look at it this way, 
# i'll have to see if there's substantial differences between doing voom 
# all at once or indiviudally
# voom_expr <- read.table(file = "ppc_rna_voom_cpm_filtered.txt",
#             sep       = "\t",
#             header = TRUE,
#             row.names = 1,
#             check.names = FALSE)
# # There is a difference between doing the voom on the whole dataframe 
# # and doing it on a per-SNP basis
# 
# voom_expr_var <- read.table(file = "ppc_rna_var_voom_cpm_filtered.txt",
#             sep       = "\t",
#             header = TRUE,
#             row.names = 1,
#             check.names = FALSE)

# Genotypes
library(snpStats)
#genotype_dosage <- fread("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/prepare_wgs/dosage_matrix.txt")
#genotype_dosage <- fread("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/phg000904.v1.NHLBI_iPSCORE_WGS.genotype-calls-vcf.c1/iPSCORE_liftover_genotypes_dosage.txt")
# genotype_dosage <- genotype_dosage %>% 
#   mutate(SNP_ID = paste0(CHROM, "_", POS, "_",REF, "_",ALT) )

# rs1635852 <- genotype_dosage[genotype_dosage$SNP_ID == "chr7_28149792_T_C",]
# rs1635852_genos <- rs1635852 %>% 
#   select(!c("CHROM", "POS", "REF", "ALT", "V273", "SNP_ID"))
# rs1635852_genos_rna <- rs1635852_genos %>% 
#   select(any_of(meta_map_ordered$WGS_UUID))
# class(rs1635852_genos_rna) <- "numeric"
# 
# rs849133 <- genotype_dosage[genotype_dosage$SNP_ID == "chr7_28152661_C_T",]
# rs849133_genos <- rs849133 %>% 
#   select(!c("CHROM", "POS", "REF", "ALT", "V273", "SNP_ID"))
# rs849133_genos_rna <- rs849133_genos %>% 
#   select(any_of(meta_map_ordered$WGS_UUID))
# class(rs849133_genos_rna) <- "numeric"
 
# rs1635852 <- genotype_dosage[genotype_dosage$SNP_ID == "chr7_28149792_T_C",]

# Set variant IDs as row names (using CHR:POS:REF:ALT)
# rownames(genotype_dosage) <- genotype_dosage$ID
# dosage <- genotype_dosage[, -(1:5)]  # drop CHROM, POS, ID, REF, ALT columns
# dosage <- t(dosage)
# class(dosage) <- "numeric"
# 
# analyzed_snps <- read_tsv(file = "analyzed_snps.txt", col_names = F)
# library(stringr)
# analyzed_snps <- str_replace(analyzed_snps$X1, "^VAR_([^_]+)_", "chr\\1_")
# 
# # extract components
# chr <- as.integer(sub("^chr([^_]+)_.*", "\\1", analyzed_snps))
# pos <- as.integer(sub("^chr[^_]+_([^_]+)_.*", "\\1", analyzed_snps))
# 
# # order
# analyzed_snps <- analyzed_snps[order(chr, pos)]
# 
# colnames(dosage) <- genotype_dosage$ID
# 
# idx <- match(analyzed_snps, colnames(dosage))
# idx <- idx[!is.na(idx)]  # drop SNPs not found
# 
# dosage_filtered <- dosage[, idx]

# # Filter: SNPs genotyped in >= 99% of individuals
# n_individuals <- nrow(dosage)
# call_rate <- colSums(!is.na(dosage)) / n_individuals
# keep_callrate <- call_rate >= 0.99
# 
# # MAF without snpStats
# af <- colMeans(dosage, na.rm = TRUE) / 2  # allele frequency of alt allele
# mafs <- pmin(af, 1 - af)                  # minor allele frequency
# keep_snps <- which(mafs > 0.05 & keep_callrate)
# dosage_filtered <- dosage[, keep_snps]

# Metadata files
# subject_meta <- read_csv("mmc2.csv")
# count_meta <- read_csv("RNA_SampleMetaData.csv")
# # 
# # # Need a map of WGS_UUID to Subject_UUID to RNA_UUID
# count_meta_ppc <- count_meta %>%
#   filter(Tissue == "PPC")
# 
# count_meta_ppc_samples <- data.frame(cbind(count_meta_ppc$WGS_UUID, count_meta_ppc$WGS_UUID))
# write_tsv(count_meta_ppc_samples, file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/ppc_samples.txt", col_names = F) 
# read_tsv(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/ppc_samples.txt")
# 
# meta_map <- count_meta_ppc %>%
#   select(RNA_UUID, WGS_UUID, Subject_UUID)
# 
# # Covariates
# # RNA Covariate File
# rna_cov_scaled <- read.table(file = "ppc_rna_covariates.txt",
#             sep       = "\t",
#             header = TRUE,
#             row.names = 1,
#             check.names = FALSE)

# lead_SNPs_PPC <- read_csv("lead_SNPs_PPC_mmc5.csv")
# lead_SNPs_PPC_eQTL <- lead_SNPs_PPC %>% 
#   filter(QTL_Type == "eQTL")
# 
# # Parse your SNP IDs into CHROM/POS format bcftools expects
# snp_df <- data.frame(
#   CHROM = sub("_.*", "", lead_SNPs_PPC_eQTL_SNPID),                        # chr10
#   POS   = as.integer(sapply(strsplit(lead_SNPs_PPC_eQTL_SNPID, "_"), `[`, 2)) # 48717378
# )
# 
# write.table(snp_df, "PPC_eQTL_lead_snps_regions.txt",
#             sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
# 
# lead_SNPs_PPC_eQTL_pos <- lead_SNPs_PPC_eQTL %>% 
#   select(Chromosome, SNP_Position)
# write_tsv(lead_SNPs_PPC_eQTL_pos, file = "lead_SNPs_PPC_eQTL_pos.tsv")
# 
# lead_SNPs_PPC_eQTL <- lead_SNPs_PPC_eQTL %>% 
#   mutate(SNP_ID = str_replace(SNP_ID, "^VAR_([^_]+)_", "chr\\1_"))
# lead_SNPs_PPC_eQTL_SNPID <- unique(lead_SNPs_PPC_eQTL$SNP_ID)
# 
# lead_SNPs_PPC_eQTL_testingPairs <- lead_SNPs_PPC_eQTL %>% 
#   select(SNP_ID, Element_ID)
# 
# # lead_SNPs_PPC_eQTL_vector <- lead_SNPs_PPC_eQTL$SNP_ID
# # lead_genes_PPC_eQTL <- lead_SNPs_PPC_eQTL$Element_ID
# 
# head(lead_SNPs_PPC_eQTL_testingPairs)

# 
# # Kinship
#kinship <- fread(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/prepare_wgs/kinship/full_kinship_v3.rel")
#kinship <- fread(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/output/kinship_v3.sXX.txt")

#kinship_id <- fread(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/prepare_wgs/kinship/to_full_kinship_v3_pruned.fam")

#kinship_id <- fread(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/prepare_wgs/kinship/to_kinship_v3_pruned.fam")
#kinship_id <- fread(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/phg000904.v1.NHLBI_iPSCORE_WGS.genotype-calls-vcf.c1/kinship/to_kinship_v2_pruned.fam")
#kinship_id <- fread(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/prepare_wgs/kinship/to_kinship_pruned_filtered.fam")
#kinship_ids <- as.vector(kinship_id$V1)
#colnames(kinship) <- kinship_ids
#rownames(kinship) <- kinship_ids

# These two structures are in the same order
# match(kinship_ids, rownames(dosage_filtered))
# 
# # Filter kinship and dosage by the WGS_UUID in the PPC tissue data
#match(kinship_ids, count_meta_ppc$WGS_UUID)
#filtered_ids <- kinship_ids[kinship_ids %in% count_meta_ppc$WGS_UUID]

# kinship <- as.matrix(kinship)
# filtered_kinship <- kinship[kinship_ids %in% rownames(genos_mat_filtered),
#                             kinship_ids %in% rownames(genos_mat_filtered)]
# 
# 
# 
# kinship_ids <- colnames(filtered_kinship)

# genos_lead <- genotype_dosage %>%
#   filter(V3 %in% lead_SNPs_PPC_eQTL_testingPairs$SNP_ID)

# genotype_dosage <- fread("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/IPSCORE_WGS_c1.hg38.snps.genotype_matrix.txt") %>%
#   rename_with(~ gsub("\\[\\d+\\]", "", .x)) %>%
#   rename_with(~ gsub("^#", "", .x)) %>%
#   rename_with(~ gsub(":GT$", "", .x))
# 
# # Convert genotypes to numeric dosage
# gt_to_dosage <- function(x) {
#   case_when(
#     x == "0/0" ~ 0L,
#     x == "0/1" ~ 1L,
#     x == "1/0" ~ 1L,  # just in case
#     x == "1/1" ~ 2L,
#     x == "./." ~ NA_integer_,
#     TRUE        ~ NA_integer_
#   )
# }
# 
# # Apply only to sample columns (not CHROM, POS, ID, REF, ALT)
# sample_cols <- colnames(genotype_dosage)[6:ncol(genotype_dosage)]
# 
# genotype_dosage <- genotype_dosage %>%
#   mutate(across(all_of(sample_cols), gt_to_dosage))
# 
# genos_mat <- as.matrix(t(genotype_dosage[,-(1:5)]))
# colnames(genos_mat) <- as.vector(genotype_dosage$ID)
# 
# genos_mat_filtered <- genos_mat[kinship_ids,]
# 
# dim(genos_mat_filtered)

# genos_lead <- genotype_dosage %>% 
#   filter(ID %in% lead_SNPs_PPC_eQTL_SNPID)

# genos_header <- read_tsv("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/prepare_wgs/dosage_matrix_header.txt")
# colnames(genos_lead) <- colnames(genos_header)
# genos_lead_mat <- as.matrix(t(genos_lead[,-(1:5)]))
# dim(genos_lead_mat)
# colnames(genos_lead_mat) <- as.character(genos_lead$ID)
# class(genos_lead_mat) = "numeric"
# colnames(genos_lead_mat) <- as.character(genos_lead$ID)
# genos_lead_SNP_IDs <- genos_lead$SNP_ID
# genos_lead_mat <- genos_lead %>% 
#   select(!c("CHROM", "POS", "REF", "ALT", "V273", "SNP_ID")) # drop CHROM, POS, ID, REF, ALT columns
# genos_lead_mat <- t(as.matrix(genos_lead_mat))
# class(genos_lead_mat) <- "numeric"
# colnames(genos_lead_mat) <- genos_lead_SNP_IDs
# 
# 
#rownames(genos_lead_mat) %in% kinship_ids
#genos_lead_mat_ordered_v2 <- genos_lead_mat[rownames(genos_lead_mat) %in% kinship_ids, ]
# save(genos_lead_mat_ordered, file = "genos_lead_mat_RNA_ordered.Rda")
#load(file = "genos_lead_mat_RNA_ordered.Rda")
#genos_lead_mat_ordered[, "chr7_28140937_T_C"]

#genos_lead_mat_ordered <- genos_lead_mat[rownames(genos_lead_mat) %in% kinship_ids,]
#save(genos_lead_mat_ordered, file = "genos_lead_SNPs_v1.Rda")
#load(file = "genos_lead_SNPs_v1.Rda")
# 
# Need to set up a pairs test dataframe such that for each gene, I'm only 
# testing SNPs that are within 100kb of it

# 
# # Using the metadata map to reorder covariates and phenotypes
# idx <- match(kinship_ids, meta_map$WGS_UUID)
# meta_map_ordered <- meta_map[idx, ]
# # 
# idx <- match(meta_map_ordered$Subject_UUID, rna_cov_scaled$Subject_UUID)
# rna_cov_ordered <- rna_cov_scaled[idx, ]
# save(rna_cov_ordered, file ="rna_cov_ordered.Rda")
# load(file ="rna_cov_ordered.Rda")
# # 
# idx <- match(meta_map_ordered$RNA_UUID, colnames(edgeR_expr))
# edgeR_expr_ordered <- edgeR_expr[, idx]
# save(edgeR_expr_ordered, file = "edgeR_expr_ordered_v2.Rda")
# load(file = "edgeR_expr_ordered_v2.Rda")

# # ====================================
# # Making the SNP-gene pair dataframe
# # ====================================
# library(GenomicRanges)
# library(dplyr)
# 
# ============================================================
# 1. BUILD SNP GRANGES
# Parse positions from SNP ID strings: chr1_918870_A_G
# ============================================================

# snp_ids <- colnames(genos_mat_filtered)
# 
# snp_gr <- data.frame(snp_id = snp_ids) |>
#   tidyr::separate(snp_id, into = c("chr", "pos", "ref", "alt"),
#                   sep = "_", remove = FALSE, extra = "drop") |>
#   mutate(pos = as.integer(pos)) |>
#   with(GRanges(
#     seqnames = chr,
#     ranges   = IRanges(start = pos, end = pos),
#     snp_id   = snp_id
#   ))
# 
# # ============================================================
# # 2. BUILD GENE GRANGES (TSS ± 1000Mb window)
# # Assumes you have a gene annotation table; adjust col names as needed
# # gene_annot should have: gene_id, chr, tss (or start, end, strand)
# # ============================================================
# 
# library(data.table)
# library(stringr)
# 
# gtf <- fread("Homo_sapiens.GRCh38.111.gtf.gz", skip = 5,
#              col.names = c("chr","source","feature","start","end",
#                            "score","strand","frame","attributes"))
# 
# # Keep only gene-level rows (one per gene, no redundant transcript/exon rows)
# gene_annot <- gtf[feature == "gene"][
#   , gene_id := str_extract(attributes, '(?<=gene_id ")[^"]+')
# ][
#   , .(gene_id, chr = paste0("chr", chr), start, end, strand)
# ]
# 
# head(gene_annot)
# 
# # Do your expr IDs have version suffixes?
# head(rownames(edgeR_expr_ordered))
# 
# # How many match as-is vs after stripping versions?
# sum(rownames(edgeR_expr_ordered) %in% gene_annot$gene_id)
# 
# # Strip version and recheck
# expr_ids_stripped <- sub("\\..*", "", rownames(edgeR_expr_ordered))
# sum(expr_ids_stripped %in% gene_annot$gene_id)
# 
# # Strip versions from expression IDs and add as a column
# gene_annot_matched <- gene_annot %>%
#   filter(gene_id %in% expr_ids_stripped) %>%
#   # Rename to make clear these are unversioned
#   dplyr::rename(gene_id_unversioned = gene_id)
# 
# # Build a lookup from unversioned -> versioned (for rejoining later)
# id_map <- data.frame(
#   gene_id_unversioned = expr_ids_stripped,
#   gene_id             = rownames(edgeR_expr_ordered)
# )
# 
# # Join versioned IDs back onto annotation
# gene_annot_matched <- gene_annot_matched |>
#   left_join(id_map, by = "gene_id_unversioned")
# 
# # Confirm
# cat("Genes with annotation:", nrow(gene_annot_matched), "\n")
# cat("Genes missing annotation:", sum(!expr_ids_stripped %in% gene_annot_matched$gene_id_unversioned), "\n")
# head(gene_annot_matched)
# 
# gene_windows <- gene_annot_matched |>
#   mutate(
#     tss       = ifelse(strand == "+", start, end),
#     win_start = pmax(1, tss - 1000000),   # 1 Mb not 100 kb
#     win_end   = tss + 1000000
#   ) |>
#   with(GRanges(
#     seqnames = chr,
#     ranges   = IRanges(start = win_start, end = win_end),
#     gene_id  = gene_id
#   ))
# 
# mhc_start <- 28510120 - 1000000
# mhc_end   <- 33480577 + 1000000
# 
# mhc_genes <- gene_annot_matched |>
#   filter(chr == "chr6",
#          end   >= mhc_start,
#          start <= mhc_end) |>
#   pull(gene_id)
# 
# cat("Genes excluded in/near MHC:", length(mhc_genes), "\n")
# 
# gene_windows <- gene_annot_matched |>
#   filter(!gene_id %in% mhc_genes) |>
#   mutate(
#     tss       = ifelse(strand == "+", start, end),
#     win_start = pmax(1, tss - 1000000),
#     win_end   = tss + 1000000
#   ) |>
#   with(GRanges(
#     seqnames = chr,
#     ranges   = IRanges(start = win_start, end = win_end),
#     gene_id  = gene_id
#   ))
# 
# # ============================================================
# # 3. FIND OVERLAPS — interval tree, scales to millions of SNPs
# # ============================================================
# 
# hits <- findOverlaps(gene_windows, snp_gr)
# 
# # ============================================================
# # 4. BUILD PAIRS DATAFRAME
# # ============================================================
# 
# testing_pairs <- data.frame(
#   Element_ID = gene_windows$gene_id[queryHits(hits)],
#   SNP_ID     = snp_gr$snp_id[subjectHits(hits)]
# )
# 
# cat("Genes tested    :", n_distinct(testing_pairs$Element_ID), "\n")
# cat("SNPs tested     :", n_distinct(testing_pairs$SNP_ID), "\n")
# cat("Total pairs     :", nrow(testing_pairs), "\n")
# cat("Avg SNPs/gene   :", round(nrow(testing_pairs) / n_distinct(testing_pairs$Element_ID)), "\n")


# ============================================================================
# SLURM entry point — pheno index comes from command line
# ============================================================================
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) < 1) stop("Usage: Rscript gemma_eqtl_analysis.R <pheno_index>")
# pheno <- as.integer(args[1])
# cat(sprintf("Processing phenotype (gene) index: %d\n", pheno))

# ---- Load pre-processed inputs ----
#load("edgeR_expr_ordered_v2.Rda")          # edgeR_expr_ordered
#load("rna_dge_filtered.Rda")                # dge_filtered
#load("GEMMA_eqtl_inputs.Rda")           # gemma_data  

# rs1635852 <- rs1635852[,-(1:5)]
# rs1635852_aligned <- rs1635852 %>% select(colnames(K))
# class(rs1635852_aligned) <- "numeric"

# Filter dge and expression matrix to only lead genes (unique)
#lead_genes   <- unique(lead_SNPs_PPC_eQTL_testingPairs$Element_ID)
#dge_lead     <- dge_filtered[, meta_map_ordered$RNA_UUID]
#expr_lead    <- edgeR_expr_ordered[rownames(edgeR_expr_ordered) %in% lead_genes, ]

# length(lead_genes)
# dim(expr_lead)
# 
# dim(genos_lead_mat_ordered)

# dim(testing_pairs)

#save(genos_mat_filtered, edgeR_expr_ordered, dge_lead, rna_cov_ordered, filtered_kinship, testing_pairs,
#     file = "allPairs_testing_inputs_050226.Rda")
load(file = "allPairs_testing_inputs_050226.Rda")

rna_cov_ordered_matrix <- rna_cov_ordered %>% 
  select(!c("Subject_UUID", "RNA_UUID"))
rna_cov_ordered_matrix <- as.matrix(rna_cov_ordered_matrix)

W <- cbind(rep(1, length(kinship_ids)), rna_cov_ordered_matrix)

#dim(W)

K <- filtered_kinship
#dim(K)

n_pairs <- nrow(testing_pairs)
head(testing_pairs)

# --- Pre-allocate results ---
def_results <- data.frame(
  SNP         = character(n_pairs),
  gene        = character(n_pairs),
  beta        = numeric(n_pairs),
  se_beta     = numeric(n_pairs),
  tau         = numeric(n_pairs),
  lambda      = numeric(n_pairs),
  F_wald      = numeric(n_pairs),
  p_wald      = numeric(n_pairs),
  stringsAsFactors = FALSE
)

# --- Pair loop ---
for (i in seq_len(n_pairs)) {
#for (i in 1:10) {
  #i <- 1
  # snp_id <- "chr7_28140937_T_C"
  #gene_id <- "ENSG00000153814.13"
  #snp_id <- "chr7_28152661_C_T"
  #snp_id <- "chr7_28149792_T_C"
  snp_id  <- testing_pairs$SNP_ID[i]
  gene_id <- testing_pairs$Element_ID[i]
  
  # Skip if SNP not in genotype matrix
  if (!snp_id %in% colnames(genos_mat_filtered)) {
    warning(sprintf("SNP %s not found in genos_lead, skipping pair %d", snp_id, i))
    next
  }
  
  if (i %% 100 == 0) cat("  Pair", i, "/", n_pairs, "\n")
  
  #causal_geno <- rs1635852_aligned 
  causal_geno  <- genos_mat_filtered[, snp_id ]          # n-vector for this SNP
  #causal_geno <- rs849133_genos_rna
  #causal_geno <- rs1635852_genos_rna
  #design_voom  <- model.matrix(~ causal_geno)
  
  # Voom on the full dge_lead (all lead genes), per-SNP design
  #y <- voom(dge_lead, design = design_voom, plot = FALSE)
  
  # Extract only the row for this gene
  expr_gene    <- edgeR_expr_ordered[gene_id, ]        # 1 × n
  #weights_gene <- y$weights[match(gene_id, rownames(dge_lead)), , drop = FALSE]  # 1 × n
  
  gemma_results <- pygemma(
    Y    = expr_gene,
    X    = causal_geno,
    W    = W,
    K    = K,
    snps = snp_id,
    verbose = 0,
    grid    = FALSE,
    nproc   = 1
  )
  #def_results_rs1635852 <- gemma_results
  
  def_results$SNP[i]     <- snp_id
  def_results$gene[i]    <- gene_id
  def_results$beta[i]    <- gemma_results$beta
  def_results$se_beta[i] <- gemma_results$se_beta
  def_results$tau[i]     <- gemma_results$tau
  def_results$lambda[i]  <- gemma_results$lambda
  def_results$F_wald[i]  <- gemma_results$F_wald
  def_results$p_wald[i]  <- gemma_results$p_wald
}
save(def_results, file = paste0("GEMMA_eqtl_output_allPairs.Rda"))
#load(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/GEMMA_eqtl_output_liftover.Rda")

# --- GEMMAmod ----

# --- Pre-allocate results ---
mod_results <- data.frame(
  SNP         = character(n_pairs),
  gene        = character(n_pairs),
  beta        = numeric(n_pairs),
  se_beta     = numeric(n_pairs),
  tau         = numeric(n_pairs),
  lambda      = numeric(n_pairs),
  F_wald      = numeric(n_pairs),
  p_wald      = numeric(n_pairs),
  stringsAsFactors = FALSE
)

# --- Pair loop ---
for (i in seq_len(n_pairs)) {
  
  snp_id  <- lead_SNPs_PPC_eQTL_testingPairs$SNP_ID[i]
  gene_id <- lead_SNPs_PPC_eQTL_testingPairs$Element_ID[i]
  
  if (!snp_id %in% colnames(genos_lead_mat_ordered)) {
    warning(sprintf("SNP %s not found in genos_lead_mat_ordered, skipping pair %d", snp_id, i))
    next
  }
  
  causal_geno <- genos_lead_mat_ordered[, snp_id]
  if (any(is.na(causal_geno))) {
    causal_geno[is.na(causal_geno)] <- mean(causal_geno, na.rm = TRUE)
  }
  
  if (i %% 100 == 0) cat("  Pair", i, "/", n_pairs, "\n")
  
  design_voom  <- model.matrix(~ causal_geno)
  
  # Voom on the full dge_lead (all lead genes), per-SNP design
  y <- voom(dge_lead, design = design_voom, plot = FALSE)
  
  # Extract only the row for this gene
  expr_gene    <- y$E[gene_id, , drop = FALSE]        # 1 × n
  weights_gene <- y$weights[match(gene_id, rownames(dge_lead)), , drop = FALSE]  # 1 × n
  
  # Shrink log-variances
  se_pheno      <- as.matrix(1 / weights_gene)
  log_s2        <- log(se_pheno)
  mu            <- mean(log_s2)
  lambda_shrink <- 0.2
  log_s2_shr    <- (1 - lambda_shrink) * log_s2 + lambda_shrink * mu
  s2_shr        <- exp(log_s2_shr)
  
  w    <- 1 / sqrt(s2_shr)
  Ww   <- diag(as.vector(w))
  
  K_whitened <- Ww %*% as.matrix(K) %*% Ww
  Y_w        <- Ww %*% t(expr_gene)
  X_w        <- Ww %*% as.matrix(causal_geno, ncol = 1)
  W_fe_w     <- Ww %*% W
  
  epsilon <- 1e-4  # slightly larger than the near-zero eigenvalue
  K_whitened_reg <- K_whitened + diag(epsilon, nrow(K_whitened))
  
  gemma_results <- pygemma(
    Y    = Y_w,
    X    = X_w,
    W    = W_fe_w,
    K    = K_whitened_reg,
    snps = snp_id,
    verbose = 0,
    grid    = FALSE,
    nproc   = 1
  )

    # Row for all markers
    mod_results$SNP[i]     <- snp_id
    mod_results$gene[i]    <- gene_id
    mod_results$beta[i]    <- gemma_results$beta
    mod_results$se_beta[i] <- gemma_results$se_beta
    mod_results$tau[i]     <- gemma_results$tau
    mod_results$lambda[i]  <- gemma_results$lambda
    mod_results$F_wald[i]  <- gemma_results$F_wald
    mod_results$p_wald[i]  <- gemma_results$p_wald
    #pval_matrix_gemma_modified[g, ] <- gemma_results$p_wald
    #if (g %% 250 == 0) cat("  Transcript", g, "/", n_genes_filt, "\n")

  }
#}
save(mod_results, file = paste0("GEMMAmod_eqtl_output_allPairs.Rda"))
#load(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/GEMMAmod_eqtl_output_liftover.Rda")

# chr7_28152661_C_T is rs849133
# chr7_28149792_T_C is rs1635852
# View(lead_SNPs_PPC[lead_SNPs_PPC$Element_Name == "JAZF1", ])
# def_results[def_results$SNP == "chr7_28140937_T_C",] # Compare the lead SNP results
# mod_results[mod_results$SNP == "chr7_28140937_T_C",]
# 
# library(ggplot2)
# library(dplyr)
# 
# # Assuming your data has columns: SNP, phenotype, model (1/2/3), pvalue
# # First, pivot wider so each model has its own p-value column
# library(tidyr)
# 
# def_results_df <- def_results %>% 
#   rename(SNP_ID = SNP,
#          Element_ID = gene) %>% 
#   mutate(log10_p_wald = -log10(p_wald), 
#          model = "GEMMA") %>% 
#   filter(!is.infinite(log10_p_wald)) %>% 
#   select(model, SNP_ID, Element_ID, log10_p_wald)
# 
# mod_results_df <- mod_results %>% 
#   rename(SNP_ID = SNP,
#          Element_ID = gene) %>% 
#   mutate(log10_p_wald = -log10(p_wald), 
#          model = "GEMMAmod") %>% 
#   filter(!is.infinite(log10_p_wald)) %>% 
#   select(model, SNP_ID, Element_ID, log10_p_wald)
# 
# lead_SNPs_PPC_eQTL_df <- lead_SNPs_PPC_eQTL %>% 
#   mutate(model = "Published",
#          log10_p_wald = -log10(P_value)) %>% 
#   select(model, SNP_ID, Element_ID, log10_p_wald)
# 
# df <- rbind(lead_SNPs_PPC_eQTL_df,
#             def_results_df)
# 
# df <- rbind(lead_SNPs_PPC_eQTL_df,
#             def_results_df,
#             mod_results_df)
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(patchwork)
# 
# # Pivot wider so each model has its own column, joined on SNP/Element pairs
# df_wide <- df %>%
#   pivot_wider(
#     id_cols = c(SNP_ID, Element_ID),
#     names_from = model,
#     values_from = log10_p_wald
#   )
# 
# # Define a consistent axis limit across all plots
# max_val <- max(df_wide$GEMMA, df_wide$GEMMAmod, df_wide$Published, na.rm = TRUE)
# 
# # Pairwise plots
# p1 <- ggplot(df_wide, aes(x = Published, y = GEMMA)) +
#   geom_point(alpha = 0.4, size = 0.8, color = "steelblue") +
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
#   coord_equal(xlim = c(0, max_val), ylim = c(0, max_val)) +
#   labs(x = "Published (-log10 p)", y = "GEMMA (-log10 p)", title = "Published vs GEMMA eQTL") +
#   theme_bw()
# p1
# 
# p2 <- ggplot(df_wide, aes(x = Published, y = GEMMAmod)) +
#   geom_point(alpha = 0.4, size = 0.8, color = "steelblue") +
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
#   coord_equal(xlim = c(0, max_val), ylim = c(0, max_val)) +
#   labs(x = "Published (-log10 p)", y = "GEMMAmod (-log10 p)", title = "Published vs GEMMAmod eQTL") +
#   theme_bw()
# 
# p3 <- ggplot(df_wide, aes(x = GEMMA, y = GEMMAmod)) +
#   geom_point(alpha = 0.4, size = 0.8, color = "steelblue") +
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
#   coord_equal(xlim = c(0, max_val), ylim = c(0, max_val)) +
#   labs(x = "GEMMA (-log10 p)", y = "GEMMAmod (-log10 p)", title = "GEMMA vs GEMMAmod eQTL") +
#   theme_bw()
# 
# p1 + p2 + p3
# 
# # Check out the characteristics of those weird SNPs in the tail, low published but not GEMMA
# df_wide_abr <- df_wide %>% 
#   filter(Published >= 5 & GEMMA <= 5)
# 
# df_wide_abr_SNPs <- df_wide_abr$SNP_ID
# 
# lead_SNPs_PPC_eQTL_abr <- lead_SNPs_PPC_eQTL %>% 
#   filter(SNP_ID %in% df_wide_abr_SNPs)
# 
# lead_SNPs_PPC_eQTL_abr_genes <- unique(lead_SNPs_PPC_eQTL_abr$Element_ID)
# 
# # Low counts?
# setwd("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/")
# 
# subject_meta <- read_csv("mmc2.csv")
# 
# counts <- read.delim("ppc_gene_counts.txt", row.names = 1)
# 
# # Need to remove leading X from any column names
# colnames(counts) <- gsub("X", "", colnames(counts))
# # Need to replace "." with "-"
# colnames(counts) <- gsub("\\.", "-", colnames(counts))
# 
# counts_abr <- counts[lead_SNPs_PPC_eQTL_abr_genes,]
# 
# sum(rowMeans(counts_abr == 0 ) > 0.20)/nrow(counts_abr) 
# 
# sum(rowMeans(counts_abr) < 10) / nrow(counts_abr)
# 
# median_raw     = apply(counts_abr, 1, median)
# 
# ## ============================================================
# ## eQTL Expression Quality Control for Discordant SNP-Gene Pairs
# ## ============================================================
# ## Inputs expected:
# ##   raw_counts  : matrix, rows = genes, cols = samples (integer counts)
# ##   norm_expr   : matrix, rows = genes, cols = samples (normalized, e.g. post-INT)
# ##   discordant_genes : character vector of gene IDs found in discordant pairs
# ##
# ## If you don't have a norm_expr yet, the script will compute a simple
# ## TMM-normalized + inverse-normal-transformed version from raw_counts.
# ## ============================================================
# 
# library(edgeR)       # for TMM normalization and CPM
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(patchwork)   # for combining plots
# 
# # ------------------------------------------------------------
# # 1. COMPUTE CPM & NORMALIZE (skip if you supply norm_expr)
# # ------------------------------------------------------------
# counts_tested <- counts[lead_genes, ]
# 
# dge <- DGEList(counts = counts_tested)
# dge <- calcNormFactors(dge, method = "TMM")
# cpm_matrix <- cpm(dge, log = FALSE)  # library-size + TMM normalized CPM
# 
# # Inverse-normal transformation per gene (rank-based)
# int_transform <- function(x) {
#   qnorm((rank(x, na.last = "keep", ties.method = "average") - 0.5) / sum(!is.na(x)))
# }
# norm_expr <- t(apply(cpm_matrix, 1, int_transform))
# 
# # ------------------------------------------------------------
# # 2. COMPUTE PER-GENE QC METRICS FOR ALL GENES
# # ------------------------------------------------------------
# 
# compute_gene_metrics <- function(raw, cpm, norm) {
#   data.frame(
#     gene_id        = rownames(raw),
#     mean_raw       = rowMeans(raw),
#     median_raw     = apply(raw, 1, median),
#     mean_cpm       = rowMeans(cpm),
#     log2_mean_cpm  = log2(rowMeans(cpm) + 1),
#     zero_fraction  = rowMeans(raw == 0),
#     cv_norm        = apply(norm, 1, sd) / (abs(apply(norm, 1, mean)) + 1e-6),
#     expr_range     = apply(norm, 1, max) - apply(norm, 1, min),
#     # skewness of normalized expression (should be ~0 after INT)
#     skewness       = apply(norm, 1, function(x) {
#       n <- length(x); m <- mean(x); s <- sd(x)
#       (sum((x - m)^3) / n) / s^3
#     }),
#     stringsAsFactors = FALSE
#   )
# }
# 
# metrics <- compute_gene_metrics(counts_tested, cpm_matrix, norm_expr)
# 
# # ------------------------------------------------------------
# # 3. FLAG LOW-QUALITY GENES
# # ------------------------------------------------------------
# # Thresholds follow GTEx/common eQTL practice — adjust as needed
# 
# metrics <- metrics |>
#   mutate(
#     flag_low_count   = mean_cpm < 1,
#     flag_high_zeros  = zero_fraction > 0.20,
#     flag_low_cv      = cv_norm < 0.05,     # barely varies across individuals
#     flag_high_skew   = abs(skewness) > 2,  # distribution still badly skewed post-INT
#     n_flags          = flag_low_count + flag_high_zeros + flag_low_cv + flag_high_skew,
#     low_quality      = n_flags >= 1,
#     group            = ifelse(gene_id %in% lead_SNPs_PPC_eQTL_abr_genes, "Discordant", "Concordant")
#   )
# 
# cat("=== Quality Flag Summary ===\n")
# metrics |>
#   group_by(group) |>
#   summarise(
#     n_genes          = n(),
#     pct_low_count    = mean(flag_low_count)  * 100,
#     pct_high_zeros   = mean(flag_high_zeros) * 100,
#     pct_low_cv       = mean(flag_low_cv)     * 100,
#     pct_high_skew    = mean(flag_high_skew)  * 100,
#     pct_any_flag     = mean(low_quality)     * 100,
#     .groups = "drop"
#   ) |>
#   print()
# 
# # Fisher's test: are flags enriched in discordant genes?
# cat("\n=== Enrichment Tests (Fisher's exact) ===\n")
# flag_cols <- c("flag_low_count", "flag_high_zeros", "flag_low_cv", "flag_high_skew", "low_quality")
# for (fl in flag_cols) {
#   tbl <- table(metrics$group, metrics[[fl]])
#   
#   # Skip if no variation in the flag (all TRUE or all FALSE)
#   if (nrow(tbl) < 2 || ncol(tbl) < 2) {
#     cat(sprintf("%-20s  skipped (no variation)\n", fl))
#     next
#   }
#   
#   ft <- fisher.test(tbl)
#   cat(sprintf("%-20s  OR = %.2f  p = %.4f\n", fl, ft$estimate, ft$p.value))
# }
# 
# # Maybe I'm comparing the wrong output, that they did something to their raw GWAS output
# # What was their discovery eQTL stats like for these SNP-gene pairs?
# PPC_discovery_eQTL_stats <- fread("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/PPC_discovery_eqtls_stats.txt")
# PPC_discovery_eQTL_stats_abr <- PPC_discovery_eQTL_stats %>% 
#   mutate(SNP_ID = str_replace(VariantID, "^VAR_([^_]+)_", "chr\\1_"))
# PPC_discovery_eQTL_stats_abr <-  PPC_discovery_eQTL_stats_abr %>% 
#   #dplyr::rename(Element_ID = ElementID) %>% 
#     right_join(lead_SNPs_PPC_eQTL_testingPairs, by = c("SNP_ID", "Element_ID"))
# head(PPC_discovery_eQTL_stats_abr)  
# 
# # The two dataframes match in their test statistics, there wasn't any additional processing...
# # WHATS HAPPENING
# # Summarize SNP characterisitcs
# lead_SNPs_PPC_eQTL_abr %>% 
#   summarise(med_MAF = median(MAF),
#             mean_MAF = mean(MAF))
# hist(lead_SNPs_PPC_eQTL_abr$MAF, xlab = "MAF", main = "MAFs of SNPs with Published << GEMMA pvals")
# # The top 10 most significant SNPs have MAF >= 0.10
# 
# # Testing their post-processing
# # eigenMT takes:
# #   --QTL      : your GEMMA output (SNP, gene, p-value columns)
# head(def_results)
# #   --GEN      : genotype matrix
# genos_lead_mat_ordered[1:5,1:5]
# #   --GENPOS   : SNP positions
# head(lead_SNPs_PPC_eQTL$SNP_ID)
# #   --PHEPOS   : gene positions (for windowing)
# head(lead_SNPs_PPC_eQTL[,4:6])
# 
# def_qtl_out <- def_results |>
#   select(
#     gene = gene,
#     snp  = SNP,
#     pvalue = p_wald
#   )
# 
# write_tsv(def_qtl_out, "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/gemma_def_qtl.txt")
# 
# # 1. Phenotypes: genes x samples
# expr_out <- as.data.frame(expr_lead)
# colnames(expr_out) <- kinship_ids
# expr_out <- tibble::rownames_to_column(expr_out, "gene_id")
# write_tsv(expr_out, "phenotypes.tsv")
# cat("phenotypes.tsv:", nrow(expr_out), "genes x", ncol(expr_out) - 1, "samples\n")
# 
# # 2. Genotypes: SNPs x samples (transpose from samples x SNPs)
# geno_out <- as.data.frame(t(genos_lead_mat_ordered))
# geno_out <- tibble::rownames_to_column(geno_out, "snp_id")
# write_tsv(geno_out, "genotypes.tsv")
# cat("genotypes.tsv:", nrow(geno_out), "SNPs x", ncol(geno_out) - 1, "samples\n")
# 
# # 3. Covariates: samples x covariates
# cov_out <- as.data.frame(W)
# colnames(cov_out) <- paste0("cov", seq_len(ncol(W)))   # name W cols first
# cov_out <- tibble::rownames_to_column(cov_out, "sample_id")  # then prepend sample_id
# cov_out$sample_id <- kinship_ids                        # overwrite with kinship_ids
# write_tsv(cov_out, "covariates.tsv")
# 
# # Should show kinship UUIDs in col 1, numeric values in remaining cols
# head(cov_out[, 1:3])
# 
# # 4. Kinship: samples x samples
# k_out <- as.data.frame(K)
# k_out <- tibble::rownames_to_column(k_out, "sample_id")
# write_tsv(k_out, "kinship.tsv")
# cat("kinship.tsv:", nrow(k_out), "x", ncol(k_out) - 1, "\n")
# 
# # 5. Pairs to test
# write_tsv(lead_SNPs_PPC_eQTL_testingPairs, "pairs.tsv")
# cat("pairs.tsv:", nrow(lead_SNPs_PPC_eQTL_testingPairs), "SNP-gene pairs\n")
# 
# limix_results <- read_tsv("limix_results.tsv")
# limix_results_filt <- limix_results %>% 
#   filter(status == "ok") %>% 
#   mutate(log10_p_lrt = -log10(p_lrt))
# limix_comp <- limix_results_filt %>% 
#   left_join(mod_results_df, by = c("SNP_ID", "Element_ID"))
# 
# plot(limix_comp$log10_p_lrt, limix_comp$log10_p_wald,
#      main = "Uncorrected -log10 p-value comparison: limix vs GEMMAmod",
#      xlab = "Limix LRT -log10 p-value",
#      ylab = "GEMMAmod Wald -log10 p-value")
# abline(a=0, b=1)
