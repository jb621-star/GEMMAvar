############################################################################
### POWER SIMULATIONS for eQTL analysis
### COMPARING: GEMMA (default) vs GEMMA (modified with voom weights)
############################################################################
.libPaths(c("/hpc/group/baughlab/jb621/R/RlibCL", .libPaths()))
library(limma)
library(edgeR)
library(tidyverse)
#library(Matrix)
library(data.table)
library(parallel)

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

############################################################################
### SIMULATION PARAMETERS
############################################################################

# Number of simulation runs
nsim <- 100

# Number of genes to simulate
ngenes <- 1000

# Number of eQTL genes (genes with true genotype effect)
n_eqtl <- 50  # 5% of genes have eQTLs
#n_eqtl <- 0 # null case 

# Number of genetic markers to simulate
n_markers <- 200

# Alpha level for FDR
alpha_fdr <- 0.05

# Load realistic abundance distribution
load(url("http://bioinf.wehi.edu.au/voom/qAbundanceDist.RData"))

# Storage for results
power_results <- data.frame(
  method = character(),
  fc = numeric(),
  sim = numeric(),
  true_positives = numeric(),
  false_positives = numeric(),
  false_negatives = numeric(),
  true_negatives = numeric(),
  power = numeric(),
  fdr_empirical = numeric(),
  stringsAsFactors = FALSE
)

################### 
# Bringing in the PLINK-formatted simulated genotypes
library(snpStats)
gbr_gen <- read.plink("/work/jb621/GEMMA_2_pickrell_sims/ped_sim/gbr_simulated_maf01")
# Convert to numeric (0, 1, 2 counts of allele 2 / minor allele)

gbr_sex <- read_tsv(file = "/work/jb621/GEMMA_2_pickrell_sims/ped_sim/gbr_samples_sex.txt", col_names = F)
# 0 = Female, 1 = Male
sex <- if_else(gbr_sex$X2 == "M", 1,0)

# Get column summaries which includes MAF
col_sum <- col.summary(gbr_gen$genotypes)

# Extract MAF specifically
mafs <- col_sum$MAF

# Keep only SNPs with MAF > 0.1
keep_snps <- which(mafs > 0.1)

# Filter the genotype matrix
geno_filtered <- gbr_gen$genotypes[, keep_snps]
map_filtered  <- gbr_gen$map[keep_snps, ]

geno_numeric <- as(geno_filtered, "numeric")

############################################################################
### RUN POWER SIMULATIONS
############################################################################

set.seed(2024)

# Randomly select 200 SNPs from the 942
# First need to prune them for MAF

tested_snps <- sample(1:ncol(geno_numeric), n_markers)
tested_snps <- sort(tested_snps)
genotype_matrix <- geno_numeric[, tested_snps]

# Olivia et al 2020
# "We identified widespread sex-biased gene expression in all tissues, 
# with 37% of genes exhibiting sex bias in at least one tissue, but with overall 
# small (median FC = 1.04) sex effects."
n_sex_biased <- 370
sex_bias_genes <- sample(1:ngenes, n_sex_biased)

kinship_matrix <- fread(file = "/work/jb621/GEMMA_2_pickrell_sims/ped_sim/output/gbr_200snps.gemmaGRM.sXX.txt")
kinship_matrix <- as.matrix(kinship_matrix)

nlibs <- ncol(kinship_matrix)

# FC parameter
fc <- 3.25 # change this parameter depending on the FC needed to run
# In our experiments we ran 1.75, 2, 2.25, 2.5, 2.75, 3, and 3.25
# fc <- 1 # null case
cat("\n==============================================\n")
cat("Testing Fold Change =", fc, "\n")
cat("==============================================\n")

for (sim in 1:nsim) {
  #sim <- 1 # for single-simulation testing
  cat("  Simulation", sim, "/", nsim, "\n")
  
  #------------------------------------------------------------------------
  # GENERATE BASELINE GENE ABUNDANCES
  #------------------------------------------------------------------------
  
  # These two, along with the BCV, are probably the most important dials of noise
  simulate_libsizes <- function(n, mean_depth = 2e7) {
    depths <- rgamma(n, shape = 1.2, scale = mean_depth / 1.2)
    outlier_idx <- sample(n, max(3, round(n * 0.2)))
    depths[outlier_idx] <- runif(length(outlier_idx), 5e4, 3e5)  # 50k-300k reads
    depths
  }
  
  expected.lib.size <- simulate_libsizes(n = nlibs)
  #hist(expected.lib.size)
  
  # Use realistic abundance distribution
  baselineprop <- qAbundanceDist((1:ngenes)/(ngenes+1))
  baselineprop <- baselineprop^2.5  # Exaggerate differences
  baselineprop <- baselineprop / sum(baselineprop)
  
  
  #------------------------------------------------------------------------
  # ASSIGN eQTL STATUS AND CAUSAL MARKERS
  #------------------------------------------------------------------------
  
  # Randomly select which genes have eQTL effects
  eqtl_genes <- sample(1:ngenes, n_eqtl)
  
  # Assign each eQTL gene to a random marker
  causal_marker <- sample(1:n_markers, n_eqtl, replace = FALSE)
  
  # Store eQTL information
  eqtl_info <- data.frame(
    gene_id = eqtl_genes,
    marker_id = causal_marker
  )
  
  # Direction of effect
  eqtl_direction <- c(rep(1, 46), rep(-1,45))
  
  #------------------------------------------------------------------------
  # GENERATE EXPECTED COUNTS WITH GENOTYPE EFFECTS (MODIFIED APPROACH)
  #------------------------------------------------------------------------
  
  # Create genotype-specific baseline proportions for each individual
  # Apply FC to proportions, not to mu0
  baselineprop_matrix <- matrix(baselineprop, ngenes, nlibs)
  
  # Apply genotype effects to baseline proportions
  for (i in 1:n_eqtl) {
    gene_idx <- eqtl_genes[i]
    marker_idx <- causal_marker[i]
    direction <- eqtl_direction[i]
    
    # Get genotype for this marker (0, 1, 2)
    geno <- genotype_matrix[, marker_idx]
    
    # Apply additive effect to baseline proportions for each individual
    for (j in 1:nlibs) {
      if (direction > 0) {
        fc_multiplier <- fc^geno[j]  # 1, fc, fc^2 for genotypes 0, 1, 2
      } else {
        fc_multiplier <- fc^(2 - geno[j])  # fc^2, fc, 1 for genotypes 0, 1, 2
      }
      baselineprop_matrix[gene_idx, j] <- baselineprop_matrix[gene_idx, j] * fc_multiplier
    }
  }
  
  # Apply the effect of sex inside loop after creating baselineprop_matrix,
  # before eQTL effects
  for (g in sex_bias_genes) {
    sex_fc    <- 1.2
    direction <- sample(c(1, -1), 1)
    for (j in 1:nlibs) {
      multiplier <- ifelse(
        (direction > 0 & sex[j] == 1) | (direction < 0 & sex[j] == 0),
        sex_fc, 1
      )
      baselineprop_matrix[g, j] <- baselineprop_matrix[g, j] * multiplier
    }
  }
  
  # Now calculate mu0 from the modified proportions
  mu0 <- baselineprop_matrix * matrix(expected.lib.size, ngenes, nlibs, byrow = TRUE)
  
  #------------------------------------------------------------------------
  # ADD BIOLOGICAL VARIATION (BCV)
  #------------------------------------------------------------------------
  
  BCV0 <- 0.5 + 1/sqrt(mu0)  # Increased variation compared to 
  df.BCV <- 10  # Lower df = more variation
  BCV <- BCV0 * sqrt(df.BCV / rchisq(ngenes, df = df.BCV))
  
  if(NCOL(BCV) == 1) BCV <- matrix(BCV, ngenes, nlibs)
  
  shape <- 1 / BCV^2
  scale <- mu0 / shape
  mu <- matrix(rgamma(ngenes * nlibs, shape = shape, scale = scale), 
               ngenes, nlibs)
  
  #------------------------------------------------------------------------
  # ADD TECHNICAL VARIATION (POISSON)
  #------------------------------------------------------------------------
  
  counts <- matrix(rpois(ngenes * nlibs, lambda = mu), ngenes, nlibs)
  
  #========================================================================
  # PREPARE DATA FOR GEMMA
  #========================================================================
  
  d0 <- DGEList(counts = counts)
  d0 <- calcNormFactors(d0)
  
  cutoff <- 1
  keep <- which(apply(cpm(d0), 1, max) > cutoff)
  d0 <- d0[keep,]
  
  cpm_tmm <- cpm(d0, normalized.lib.sizes = TRUE)
  
  invnorm <- function(x) {
    r <- rank(x, ties.method = "average")
    qnorm((r - 0.5) / length(r))
  }
  
  expr_int <- t(apply(cpm_tmm, 1, invnorm))

  # Update eQTL info for filtered genes
  eqtl_info_filt <- eqtl_info[eqtl_info$gene_id %in% keep, ]
  # Remap gene IDs to filtered indices
  gene_id_mapping <- setNames(1:length(keep), keep)
  eqtl_info_filt$gene_id_filt <- gene_id_mapping[as.character(eqtl_info_filt$gene_id)]
  
  n_genes_filt <- length(gene_id_mapping)
  
  # Fixed effect covariates (just intercept)
  W <- matrix(1, nrow = nlibs, ncol = 1)
  W <- cbind(W, sex)
  #W <- model.matrix(~sex)
  
  # SNP names
  snp_names <- 1:n_markers
  
  #========================================================================
  # METHOD 1: GEMMA (DEFAULT - no heteroscedasticity correction)
  #========================================================================
  
  # Storage for all gene-marker p-values
  pval_matrix_gemma_default <- matrix(1, nrow = n_genes_filt, ncol = n_markers)
  
  # mm0 <- model.matrix(~ 1, data = data.frame(sample = 1:nlibs))
  # y <- voom(d0, mm0, plot = FALSE)
  
  # The reverse in this case: test marker against all genes
  for (m in 1:n_markers) {
    #m <- 1
    if (m %% 50 == 0) cat("  SNP", m, "/", n_markers, "\n")
    
    # causal_geno <- genotype_matrix[, m]
    # design_voom <- model.matrix(~ causal_geno)
    # 
    # y <- voom(d0, design = design_voom, plot = F)
    
    for (g in 1:n_genes_filt) {
      #g <- 1
      # Get voom-transformed expression with appropriate design
      # expr_gene <- y$E[g, , drop = FALSE]
      
      # TMM-based
      expr_gene <- expr_int[g, , drop = FALSE]
      
      # Prepare data for GEMMA
      data <- prepare_gemma_data(
        pheno = expr_gene,
        geno = genotype_matrix[, m],
        #geno = genotype_matrix,
        covars = W,
        kinship = kinship_matrix
      )
      
      # Run GEMMA (default)
      gemma_results <- pygemma(
        Y = data$Y,
        X = data$X,
        W = data$W,
        K = data$K,
        snps = snp_names[m],
        verbose = 0,
        grid = FALSE,
        nproc = 4
      )
      
      # Extract p-values for all markers
      pval_matrix_gemma_default[g, m] <- gemma_results$p_wald
      
    }
  }
  
  # Apply FDR correction across all tests
  
  adj_pvals_default <- p.adjust(as.vector(pval_matrix_gemma_default), method = "BH")
  #adj_pvals_default <- as.vector(pval_matrix_gemma_default)
  adj_pval_matrix_default <- matrix(adj_pvals_default, nrow = n_genes_filt, ncol = n_markers)
  
  # # A gene is called significant if its best marker passes FDR threshold
  sig_genes_default <- apply(adj_pval_matrix_default, 1, min) < alpha_fdr
  
  # For each gene, find the most significant marker
  best_marker_default <- apply(adj_pval_matrix_default, 1, which.min)
  best_pval_default <- apply(adj_pval_matrix_default, 1, min)
  
  #========================================================================
  # METHOD 2: GEMMA (MODIFIED - with voom-like weights)
  #========================================================================
  
  # Storage for all gene-marker p-values
  pval_matrix_gemma_modified <- matrix(1, nrow = n_genes_filt, ncol = n_markers)
  
  # The reverse in this case: test marker against all genes
  for (m in 1:n_markers) {
    #m <- 9 # ind 54 and 18 are outliers
    if (m %% 50 == 0) cat("  SNP", m, "/", n_markers, "\n")
    causal_geno <- genotype_matrix[, m]
    design_voom <- model.matrix(~ causal_geno)
    
    y <- voom(d0, design = design_voom, plot = F)
    
    for (g in 1:n_genes_filt) {
      # Get voom-transformed expression with appropriate design
      expr_gene <- y$E[g, , drop = FALSE]
      weights_gene <- y$weights[g, , drop = FALSE]
      
      # Calculate variance estimates from voom weights
      se_pheno <- as.matrix(1/weights_gene)
      log_s2 <- log(se_pheno)
      mu <- mean(log_s2)
      
      lambda_shrink <- 0
      log_s2_shr <- (1 - lambda_shrink) * log_s2 + lambda_shrink * mu
      s2_shr <- exp(log_s2_shr)
      
      w <- 1 / sqrt(s2_shr)
      Ww <- diag(as.vector(w))
      
      # Whiten
      K_whitened <- Ww %*% as.matrix(kinship_matrix) %*% Ww
      
      Y_w <- Ww %*% t(expr_gene)
      X_w <- Ww %*% genotype_matrix             # X is n × p
      W_fe_w <- Ww %*% W           # fixed effects
      
      # Prepare data
      data <- prepare_gemma_data(
        pheno = Y_w,
        geno = X_w[, m],
        covars = W_fe_w,
        kinship = K_whitened
      )
      
      # Run analysis
      gemma_results <- pygemma(
        Y = data$Y,
        X = data$X,
        W = W_fe_w,
        K = data$K,
        snps = snp_names[m],
        verbose = 0,
        grid = FALSE,  # Use fast optimization
        nproc = 4      # Use 4 cores
      )
      
      # Extract p-values for all markers
      pval_matrix_gemma_modified[g, m] <- gemma_results$p_wald

    }
  }
  
  # For each gene, find the most significant marker
  # Apply FDR correction across all tests
  adj_pvals_modified <- p.adjust(as.vector(pval_matrix_gemma_modified), method = "BH")
  #adj_pvals_modified <- as.vector(pval_matrix_gemma_modified)
  adj_pval_matrix_modified <- matrix(adj_pvals_modified, nrow = n_genes_filt, ncol = n_markers)
  
  # A gene is called significant if its best marker passes FDR threshold
  sig_genes_modified <- apply(adj_pval_matrix_modified, 1, min) < alpha_fdr
  
  best_marker_modified <- apply(adj_pval_matrix_modified, 1, which.min)
  best_pval_modified <- apply(adj_pval_matrix_modified, 1, min)
  
  #========================================================================
  # EVALUATE PERFORMANCE FOR BOTH METHODS
  #========================================================================
  
  # Create truth vectors
  is_eqtl <- rep(FALSE, n_genes_filt)
  true_marker <- rep(NA, n_genes_filt)
  
  for (i in 1:nrow(eqtl_info_filt)) {
    gene_idx <- eqtl_info_filt$gene_id_filt[i]
    is_eqtl[gene_idx] <- TRUE
    true_marker[gene_idx] <- eqtl_info_filt$marker_id[i]
  }
  
  # ------- GEMMA DEFAULT -------
  TP_default <- 0
  FP_default <- 0
  FN_default <- 0
  TN_default <- 0
  
  for (g in 1:n_genes_filt) {
    if (is_eqtl[g]) {
      if (sig_genes_default[g]) {
        if (best_marker_default[g] == true_marker[g]) {
          TP_default <- TP_default + 1  # Correct detection
        } else {
          FP_default <- FP_default + 1  # Wrong marker
        }
      } else {
        FN_default <- FN_default + 1  # Missed the eQTL
      }
    } else {
      if (sig_genes_default[g]) {
        FP_default <- FP_default + 1  # False positive
      } else {
        TN_default <- TN_default + 1  # Correct null
      }
    }
  }
  
  power_default <- ifelse(sum(is_eqtl) > 0, TP_default / sum(is_eqtl), 0)
  fdr_emp_default <- ifelse(TP_default + FP_default > 0,
                            FP_default / (TP_default + FP_default), 0)
  
  # ------- GEMMA DEFAULT - NULL CASE -------
  # FP_default <- sum(sig_genes_default)   # All detections are FPs
  # TN_default <- sum(!sig_genes_default)  # All non-detections are TNs
  # TP_default <- 0
  # FN_default <- 0
  # 
  # power_default  <- NA_real_             # Undefined - no true positives exist
  # fdr_emp_default <- ifelse(FP_default > 0, 1.0, 0.0)  # All discoveries are false
  # fpr_default <- FP_default / n_genes_filt              # Type I error rate
  # 
  # Store results - GEMMA DEFAULT
  power_results <- rbind(power_results, data.frame(
    method = "GEMMA_default",
    fc = fc,
    sim = sim,
    true_positives = TP_default,
    false_positives = FP_default,
    false_negatives = FN_default,
    true_negatives = TN_default,
    power = power_default,
    fdr_empirical = fdr_emp_default,
    stringsAsFactors = FALSE
  ))
  
  # ------- GEMMA MODIFIED -------
  TP_modified <- 0
  FP_modified <- 0
  FN_modified <- 0
  TN_modified <- 0
  
  for (g in 1:n_genes_filt) {
    if (is_eqtl[g]) {
      if (sig_genes_modified[g]) {
        if (best_marker_modified[g] == true_marker[g]) {
          TP_modified <- TP_modified + 1  # Correct detection
        } else {
          FP_modified <- FP_modified + 1  # Wrong marker
        }
      } else {
        FN_modified <- FN_modified + 1  # Missed the eQTL
      }
    } else {
      if (sig_genes_modified[g]) {
        FP_modified <- FP_modified + 1  # False positive
      } else {
        TN_modified <- TN_modified + 1  # Correct null
      }
    }
  }
  
  power_modified <- ifelse(sum(is_eqtl) > 0, TP_modified / sum(is_eqtl), 0)
  fdr_emp_modified <- ifelse(TP_modified + FP_modified > 0,
                             FP_modified / (TP_modified + FP_modified), 0)
  
  # # ------- GEMMA MODIFIED - NULL CASE -------
  # # FP_default <- sum(sig_genes_default)   # All detections are FPs
  # # TN_default <- sum(!sig_genes_default)  # All non-detections are TNs
  # # TP_default <- 0
  # # FN_default <- 0
  # # 
  # # power_default  <- NA_real_             # Undefined - no true positives exist
  # # fdr_emp_default <- ifelse(FP_default > 0, 1.0, 0.0)  # All discoveries are false
  # # fpr_default <- FP_default / n_genes_filt              # Type I error rate
  # # 
  # # Store results - GEMMA MODIFIED
  power_results <- rbind(power_results, data.frame(
    method = "GEMMA_modified",
    fc = fc,
    sim = sim,
    true_positives = TP_modified,
    false_positives = FP_modified,
    false_negatives = FN_modified,
    true_negatives = TN_modified,
    power = power_modified,
    fdr_empirical = fdr_emp_modified,
    stringsAsFactors = FALSE
  ))
}
#}

############################################################################
### SUMMARIZE RESULTS
############################################################################

# Calculate mean power and FDR for each method and fold change
summary_stats <- power_results %>%
  group_by(method, fc) %>%
  summarise(
    mean_power = mean(power, na.rm = TRUE),
    sd_power = sd(power, na.rm = TRUE),
    mean_fdr = mean(fdr_empirical, na.rm = TRUE),
    sd_fdr = sd(fdr_empirical, na.rm = TRUE),
    mean_TP = mean(true_positives),
    mean_FP = mean(false_positives),
    mean_FN = mean(false_negatives),
    mean_TN = mean(true_negatives),
    .groups = "drop"
  )

print(summary_stats)

# Create comparison plots
# library(ggplot2)
# 
# # Power comparison
# p1 <- ggplot(summary_ci, aes(x = fc, y = power_pooled, color = method, group = method)) +
#   geom_line(size = 1) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = power_lower, ymax = power_upper),
#                 width = 0.02) +
#   scale_color_discrete(labels = c("GEMMA_default" = "GEMMA",
#                                   "GEMMA_modified" = "GEMMAmod")) +
#   labs(title = "True Detection Rate Comparison: GEMMA vs GEMMAmod",
#        subtitle = "Modified version uses voom-derived variances for heteroscedasticity",
#        x = "Fold Change of Causal SNP",
#        y = "TDR (Correct eQTL Detection)",
#        color = "Method") +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   theme(axis.text.x=element_text(size=14),
#         axis.text.y=element_text(size=14),
#         plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
#         plot.subtitle = element_text(hjust = 0.5, size = 14),
#         legend.title = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         axis.title.x = element_text(size = 18),             # Increase x-axis label size
#         axis.title.y = element_text(size = 18))             # Increase y-axis label size
# p1
# 
# # FDR comparison
# p2 <- ggplot(summary_ci, aes(x = fc, y = fdr_pooled, color = method, group = method)) +
#   geom_line(size = 1) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = fdr_lower, ymax = fdr_upper),
#                 width = 0.02) +
#   scale_color_discrete(labels = c("GEMMA_default" = "GEMMA",
#                                   "GEMMA_modified" = "GEMMAmod")) +
#   geom_hline(yintercept = alpha_fdr, linetype = "dashed", color = "red") +
#   labs(title = "FDR Control: GEMMA vs GEMMAmod",
#        subtitle = "Red line shows nominal FDR threshold (0.05)",
#        x = "Fold Change",
#        y = "Empirical FDR",
#        color = "Method") +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   theme(axis.text.x=element_text(size=14),
#         axis.text.y=element_text(size=14),
#         plot.title = element_text(hjust = 0.5, size = 22),  # Center the title and set its size
#         plot.subtitle = element_text(hjust = 0.5, size = 14),
#         legend.title = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         axis.title.x = element_text(size = 18),             # Increase x-axis label size
#         axis.title.y = element_text(size = 18))             # Increase y-axis label size
# p2
# 
# # Save plots
# ggsave("power_comparison_gemma.png", p1, width = 10, height = 6)
# ggsave("fdr_comparison_gemma.png", p2, width = 10, height = 6)
#save(pval_matrix_gemma_default, pval_matrix_gemma_modified, file = "/work/jb621/GEMMA_2_pickrell_sims/eqtl_gemma_fc15_power_simulation_results_022526.RData")

# Save results
save(power_results, summary_stats, kinship_matrix,
     file = "/work/jb621/GEMMA_2_pickrell_sims/eqtl_gemma_voomVedgeR_NOSHRINKAGE_power_simulation_fc325_results_041326.RData")

