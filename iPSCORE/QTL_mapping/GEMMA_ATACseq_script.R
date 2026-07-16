.libPaths(c("/home/jb621/R/x86_64-pc-linux-gnu-library/4.2","/usr/local/lib/R/site-library","/usr/local/lib/R/library",.libPaths()) ) 
library(limma)
library(edgeR)
library(tidyverse)
#library(Matrix)
library(data.table)
library(parallel)

setwd("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/")

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
# Normalized + filtered + INT-transformed matrix (for eQTL input)
edgeR_expr <- read.table(file = "ppc_atac_tmm_filtered.txt",
            sep       = "\t",
            header = TRUE,
            row.names = 1,
            check.names = FALSE)   # col.names = NA writes row names correctly

# VOOM CPM (post-filter) 
# it might be incorrect to look at it this way, 
# i'll have to see if there's substantial differences between doing voom 
# all at once or indiviudally
voom_expr <- read.table(file = "ppc_atac_voom_cpm_filtered.txt",
            sep       = "\t",
            header = TRUE,
            row.names = 1,
            check.names = FALSE)
# There is a difference between doing the voom on the whole dataframe
# and doing it on a per-SNP basis

voom_expr_var <- read.table(file = "ppc_atac_var_voom_cpm_filtered.txt",
            sep       = "\t",
            header = TRUE,
            row.names = 1,
            check.names = FALSE)

# Genotypes
library(snpStats)
genotype_dosage <- fread("dosage_matrix_final.txt")
rs1635852 <- genotype_dosage[genotype_dosage$ID == "chr7_28149792_T_C",]

# Set variant IDs as row names (using CHR:POS:REF:ALT)
# rownames(genotype_dosage) <- genotype_dosage$ID
dosage <- genotype_dosage[, -(1:5)]  # drop CHROM, POS, ID, REF, ALT columns
dosage <- t(dosage)
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
colnames(dosage) <- genotype_dosage$ID
WGS_IDs <- colnames(K)
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
subject_meta <- read_csv("mmc2.csv")
count_meta <- read_csv("ATAC_SampleMetaData.csv")

# Need a map of WGS_UUID to Subject_UUID to RNA_UUID
count_meta_ppc <- count_meta %>%
  filter(Tissue == "PPC")
meta_map <- count_meta_ppc %>%
  select(ATAC_UUID, WGS_UUID, Subject_UUID)

# # Covariates
# # RNA Covariate File
atac_cov_scaled <- read.table(file = "ppc_ATAC_covariates.txt",
            sep       = "\t",
            header = TRUE,
            row.names = 1,
            check.names = FALSE)
# 
# # Kinship
# kinship <- fread(file = "kinship_full.sXX.txt")
# 
# kinship_id <- fread("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/prepare_wgs/kinship/to_kinship_pruned.fam")
# kinship_ids <- as.vector(kinship_id$V1)
# colnames(kinship) <- kinship_ids
# rownames(kinship) <- kinship_ids
# 
# # These two structures are in the same order
# #match(kinship_ids, rownames(dosage_filtered))
# 
# # Filter kinship and dosage by the WGS_UUID in the PPC tissue data
# #match(kinship_ids, count_meta_ppc$WGS_UUID)
# filtered_ids <- kinship_ids[kinship_ids %in% count_meta_ppc$WGS_UUID]
# 
# kinship <- as.matrix(kinship)
# filtered_kinship <- kinship[kinship_ids %in% count_meta_ppc$WGS_UUID, 
#                             kinship_ids %in% count_meta_ppc$WGS_UUID]
# 
# dosage_filtered_2 <- dosage_filtered[kinship_ids %in% count_meta_ppc$WGS_UUID,]
# dosage_filtered_2 <- matrix(as.numeric(dosage_filtered_2), ncol = 5323793)
# colnames(dosage_filtered_2) <- colnames(dosage_filtered)
# 
# # Using the metadata map to reorder covariates and phenotypes
idx <- match(colnames(K), meta_map$WGS_UUID)
meta_map_ordered <- meta_map[idx, ]
# 
idx <- match(meta_map_ordered$Subject_UUID, atac_cov_scaled$Subject_UUID)
atac_cov_ordered <- atac_cov_scaled[idx, ]
W <- atac_cov_ordered %>% 
  select(!c("Subject_UUID", "ATAC_UUID"))
W <- as.matrix(W)

idx <- match(meta_map_ordered$ATAC_UUID, colnames(edgeR_expr))
edgeR_expr_ordered <- edgeR_expr[, idx]
save(edgeR_expr_ordered, file = "edgeR_atac_ordered.Rda")
#load(file = "edgeR_atac_ordered.Rda")

# ============================================================================
# SLURM entry point — pheno index comes from command line
# ============================================================================

# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) < 1) stop("Usage: Rscript gemma_eqtl_analysis.R <pheno_index>")
# pheno <- as.integer(args[1])
# cat(sprintf("Processing phenotype (gene) index: %d\n", pheno))

# ---- Load pre-processed inputs ----
load("edgeR_atac_ordered.Rda")          # edgeR_expr_ordered
load("atac_dge_filtered.Rda")                # dge_filtered

lead_SNPs_PPC <- read_csv("lead_SNPs_PPC_mmc5.csv")
lead_SNPs_PPC_caQTL <- lead_SNPs_PPC %>% 
  filter(QTL_Type == "caQTL")

lead_SNPs_PPC_caQTL_pos <- lead_SNPs_PPC_caQTL %>% 
  select(Chromosome, SNP_Position)
write_tsv(lead_SNPs_PPC_caQTL_pos, file = "lead_SNPs_PPC_caQTL_pos.tsv")

lead_SNPs_PPC_caQTL <- lead_SNPs_PPC_caQTL %>% 
  mutate(SNP_ID = str_replace(SNP_ID, "^VAR_([^_]+)_", "chr\\1_"))

lead_SNPs_PPC_caQTL_testingPairs <- lead_SNPs_PPC_caQTL %>% 
  select(SNP_ID, Element_ID)

# lead_SNPs_PPC_caQTL_vector <- lead_SNPs_PPC_caQTL$SNP_ID
# lead_genes_PPC_caQTL <- lead_SNPs_PPC_caQTL$Element_ID

head(lead_SNPs_PPC_caQTL_testingPairs)
edgeR_expr_ordered[1:5,1:5]
head(dge_aligned)

genos_lead <- dosage[ , colnames(dosage) %in% lead_SNPs_PPC_caQTL_testingPairs$SNP_ID]
genos_lead <- genos_lead[rownames(genos_lead) %in% colnames(K),]
class(genos_lead) <- "numeric"
save(genos_lead, file = "genos_atac_lead_SNPs.Rda")

rs1635852 <- rs1635852[,-(1:5)]
rs1635852_aligned <- rs1635852 %>% select(colnames(K))
class(rs1635852_aligned) <- "numeric"

# Filter dge and expression matrix to only lead genes (unique)
lead_peaks   <- unique(lead_SNPs_PPC_caQTL_testingPairs$Element_ID)
dge_lead     <- dge_aligned[rownames(dge_aligned) %in% lead_peaks, ]
expr_lead    <- edgeR_expr_ordered[rownames(edgeR_expr_ordered) %in% lead_peaks, ]

K <- kinship

n_pairs <- nrow(lead_SNPs_PPC_caQTL_testingPairs)

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
  #i <- 1
  #snp_id <- "chr7_28149792_T_C"
  #gene_id <- "ENSG00000153814.13"
  snp_id  <- lead_SNPs_PPC_caQTL_testingPairs$SNP_ID[i]
  gene_id <- lead_SNPs_PPC_caQTL_testingPairs$Element_ID[i]
  
  # Skip if SNP not in genotype matrix
  if (!snp_id %in% colnames(genos_lead)) {
    warning(sprintf("SNP %s not found in genos_lead, skipping pair %d", snp_id, i))
    next
  }
  
  if (i %% 100 == 0) cat("  Pair", i, "/", n_pairs, "\n")
  
  #causal_geno <- rs1635852_aligned 
  causal_geno  <- genos_lead[, snp_id ]          # n-vector for this SNP
  #design_voom  <- model.matrix(~ causal_geno)
  
  # Voom on the full dge_lead (all lead peaks), per-SNP design
  #y <- voom(dge_lead, design = design_voom, plot = FALSE)
  
  # Extract only the row for this gene
  expr_gene    <- expr_lead[gene_id, ]        # 1 × n
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
save(def_results, file = paste0("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/GEMMA_caQTL_output.Rda"))
load(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/GEMMA_caQTL_output.Rda")

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
  
  snp_id  <- lead_SNPs_PPC_caQTL_testingPairs$SNP_ID[i]
  gene_id <- lead_SNPs_PPC_caQTL_testingPairs$Element_ID[i]
  
  if (!snp_id %in% colnames(genos_lead)) {
    warning(sprintf("SNP %s not found in genos_lead, skipping pair %d", snp_id, i))
    next
  }
  
  causal_geno <- as.matrix(genos_lead[, snp_id])
  if (any(is.na(causal_geno))) {
    causal_geno[is.na(causal_geno)] <- mean(causal_geno, na.rm = TRUE)
  }
  
  if (i %% 100 == 0) cat("  Pair", i, "/", n_pairs, "\n")
  
  design_voom  <- model.matrix(~ causal_geno)
  
  # Voom on the full dge_lead (all lead peaks), per-SNP design
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
  
  gemma_results <- pygemma(
    Y    = Y_w,
    X    = X_w,
    W    = W_fe_w,
    K    = K_whitened,
    snps = snp_id,
    verbose = 0,
    grid    = FALSE,
    nproc   = 1
  )

    # Row for all markers
    mod_results$SNP[i]     <- gemma_results$SNP
    mod_results$beta[i]    <- gemma_results$beta
    mod_results$se_beta[i] <- gemma_results$se_beta
    mod_results$tau[i]     <- gemma_results$tau
    mod_results$lambda[i]  <- gemma_results$lambda
    mod_results$F_wald[i]  <- gemma_results$F_wald
    mod_results$p_wald[i]  <- gemma_results$p_wald
    #pval_matrix_gemma_modified[g, ] <- gemma_results$p_wald
    #if (g %% 250 == 0) cat("  Transcript", g, "/", n_peaks_filt, "\n")

  }
#}
save(mod_results, file = paste0("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/GEMMAmod_caQTL_output.Rda"))
load(file = "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/GEMMAmod_caQTL_output.Rda")

# rs864745 is chr7_28140937_T_C
lead_SNPs_PPC[lead_SNPs_PPC$Element_Name == "JAZF1", "P_value"]
def_results[def_results$SNP == "chr7_28140937_T_C",]
mod_results[mod_results$SNP == "chr7_28140937_T_C",]

plot(def_results$p_wald, mod_results$p_wald)
abline(a=0, b=1)

# Single test for rs1635852 on JAZF1
def_results_rs1635852
mod_results_rs1635852

# The two SNPs are in perfect LD


