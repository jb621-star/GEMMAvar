.libPaths(c("/home/jb621/R/x86_64-pc-linux-gnu-library/4.2","/usr/local/lib/R/site-library","/usr/local/lib/R/library",.libPaths()) ) 
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

setwd("/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping")

n_sims <- 1000

#fread(file = "/hpc/group/baughlab/jb621/GEMMA_2_testing/slope100/output/gemmaGRM.NORMAL.sXX.txt") -> kinship2
#load(file = "/work/jb621/simHaplo/LD_mapping/taggingSNPs_v4_geno_df.Rda")
# ============================================================================
# Power Analysis for GEMMA Methods
# ============================================================================
# Functions to calculate and compare power between Default and Modified GEMMA
# across different heritability levels
# ============================================================================

# ============================================================================
# File Reading Functions
# ============================================================================

#' Read a single phenotype file
#'
#' @param file_path Path to phenotype file
#' @return Data frame or vector of phenotype values
read_pheno_file <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(sprintf("File not found: %s", file_path))
    return(NULL)
  }
  
  tryCatch({
    # Try reading as table first (may have headers/row names)
    data <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
    
    # If multiple columns, assume first is strain ID, second is phenotype
    if (ncol(data) >= 2) {
      colnames(data)[1:2] <- c("strain", "value")
      return(data[, 1:2])
    } else if (ncol(data) == 1) {
      # Single column - need to add strain IDs
      colnames(data) <- "value"
      return(data)
    }
  }, error = function(e) {
    warning(sprintf("Error reading %s: %s", file_path, e$message))
    return(NULL)
  })
}

#' Read a single variance file
#'
#' @param file_path Path to variance file
#' @return Data frame or vector of variance values
read_se2_file <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(sprintf("File not found: %s", file_path))
    return(NULL)
  }
  
  tryCatch({
    data <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
    
    if (ncol(data) >= 2) {
      colnames(data)[1:2] <- c("strain", "se2")
      return(data[, 1:2])
    } else if (ncol(data) == 1) {
      colnames(data) <- "se2"
      return(data)
    }
  }, error = function(e) {
    warning(sprintf("Error reading %s: %s", file_path, e$message))
    return(NULL)
  })
}

#' Read and compile a single simulation
#'
#' @param sim_id Simulation ID (1 to 1000)
#' @param pheno_template Template for phenotype file path
#' @param se2_template Template for variance file path
#' @param strain_names Optional vector of strain names
#' @return Data frame with simulation data
read_single_simulation <- function(sim_id, pheno_template, se2_template, 
                                   strain_names = NULL) {
  
  # Construct file paths
  pheno_file <- sprintf(pheno_template, sim_id)
  se2_file <- sprintf(se2_template, sim_id)
  
  # Read files
  pheno_data <- read_pheno_file(pheno_file)
  se2_data <- read_se2_file(se2_file)
  
  # Check if both files were read successfully
  if (is.null(pheno_data) || is.null(se2_data)) {
    return(NULL)
  }
  
  # Add strain names if not present and provided
  if (!"strain" %in% colnames(pheno_data) && !is.null(strain_names)) {
    pheno_data$strain <- strain_names
  }
  if (!"strain" %in% colnames(se2_data) && !is.null(strain_names)) {
    se2_data$strain <- strain_names
  }
  
  # Merge phenotype and variance data
  if ("strain" %in% colnames(pheno_data) && "strain" %in% colnames(se2_data)) {
    combined <- inner_join(pheno_data, se2_data, by = "strain")
  } else {
    # If no strain column, assume same order
    combined <- cbind(pheno_data, se2_data)
    if (!"strain" %in% colnames(combined) && !is.null(strain_names)) {
      combined$strain <- strain_names
    }
  }
  
  # Add simulation metadata
  combined$sim_id <- sim_id
  combined$trait <- sprintf("Trait_%d", sim_id)
  
  return(combined)
}

# ============================================================================
# Main Compilation Function
# ============================================================================

#' Compile all simulated phenotype files
#'
#' @param n_sims Number of simulations to read
#' @param pheno_template Template for phenotype file paths
#' @param se2_template Template for variance file paths
#' @param strain_names Optional vector of strain names
#' @param n_cores Number of cores for parallel processing
#' @param sample_first_file Logical, whether to examine first file for structure
#' @return Data frame with all compiled data
compile_simulated_data <- function(n_sims, pheno_template, se2_template,
                                   strain_names = NULL, n_cores = 1,
                                   sample_first_file = TRUE) {
  
  cat(sprintf("Compiling %d simulated phenotype files...\n", n_sims))
  
  # First, examine the structure of the first file
  if (sample_first_file) {
    cat("\n=== Examining first file structure ===\n")
    first_pheno <- read_pheno_file(sprintf(pheno_template, 1))
    first_se2 <- read_se2_file(sprintf(se2_template, 1))
    
    cat("\nPhenotype file structure:\n")
    print(head(first_pheno))
    cat(sprintf("Dimensions: %d rows × %d columns\n", nrow(first_pheno), ncol(first_pheno)))
    
    cat("\nVariance file structure:\n")
    print(head(first_se2))
    cat(sprintf("Dimensions: %d rows × %d columns\n", nrow(first_se2), ncol(first_se2)))
    
    # Extract strain names if present in first file
    if (is.null(strain_names) && "strain" %in% colnames(first_pheno)) {
      strain_names <- first_pheno$strain
      cat(sprintf("\nExtracted %d strain names from first file\n", length(strain_names)))
    }
  }
  
  # Read all simulations
  cat("\n=== Reading all simulation files ===\n")
  
  if (n_cores > 1) {
    cat(sprintf("Using %d cores for parallel processing...\n", n_cores))
    
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    
    # Export necessary objects
    clusterExport(cl, c("pheno_template", "se2_template", "strain_names",
                        "read_pheno_file", "read_se2_file", 
                        "read_single_simulation"),
                  envir = environment())
    
    # Load required packages on each worker
    clusterEvalQ(cl, library(tidyverse))
    
    # Read in parallel
    if (requireNamespace("pbapply", quietly = TRUE)) {
      results_list <- pbapply::pblapply(1:n_sims, function(sim_id) {
        read_single_simulation(sim_id, pheno_template, se2_template, strain_names)
      }, cl = cl)
    } else {
      results_list <- parLapply(cl, 1:n_sims, function(sim_id) {
        read_single_simulation(sim_id, pheno_template, se2_template, strain_names)
      })
    }
  } else {
    # Sequential processing
    if (requireNamespace("pbapply", quietly = TRUE)) {
      results_list <- pbapply::pblapply(1:n_sims, function(sim_id) {
        read_single_simulation(sim_id, pheno_template, se2_template, strain_names)
      })
    } else {
      results_list <- lapply(1:n_sims, function(sim_id) {
        if (sim_id %% 100 == 0) cat(sprintf("  Progress: %d/%d\n", sim_id, n_sims))
        read_single_simulation(sim_id, pheno_template, se2_template, strain_names)
      })
    }
  }
  
  # Check for failed reads
  n_failed <- sum(sapply(results_list, is.null))
  if (n_failed > 0) {
    cat(sprintf("\nWarning: %d files failed to read\n", n_failed))
    failed_ids <- which(sapply(results_list, is.null))
    cat("Failed simulation IDs:", paste(head(failed_ids, 20), collapse = ", "))
    if (n_failed > 20) cat(sprintf(" ... and %d more", n_failed - 20))
    cat("\n")
  }
  
  # Remove NULL entries
  results_list <- results_list[!sapply(results_list, is.null)]
  
  # Combine all results
  cat("\n=== Combining all data ===\n")
  combined_data <- bind_rows(results_list)
  
  # Ensure consistent column names and structure
  # Match the structure of normal_null_data
  if (!"causal_snp" %in% colnames(combined_data)) {
    # You'll need to add this information - for now using placeholder
    combined_data$causal_snp <- NA_character_
  }
  
  if (!"method" %in% colnames(combined_data)) {
    combined_data$method <- "Simulated_Null"
  }
  
  # Reorder columns to match normal_null_data structure
  # strain, trait, causal_snp, sim_id (like replicate), value, se2, method
  column_order <- c("strain", "trait", "causal_snp", "sim_id", 
                    "value", "se2", "method")
  
  # Only select columns that exist
  column_order <- column_order[column_order %in% colnames(combined_data)]
  combined_data <- combined_data %>% select(all_of(column_order))
  
  # Print summary
  cat("\n=== Summary ===\n")
  cat(sprintf("Total observations: %s\n", format(nrow(combined_data), big.mark = ",")))
  cat(sprintf("Unique strains: %d\n", length(unique(combined_data$strain))))
  cat(sprintf("Unique traits: %d\n", length(unique(combined_data$trait))))
  cat(sprintf("Observations per trait: %d\n", 
              nrow(combined_data) / length(unique(combined_data$trait))))
  
  cat("\nData structure:\n")
  print(str(combined_data))
  
  cat("\nFirst few rows:\n")
  print(head(combined_data, 10))
  
  cat("\nSummary statistics:\n")
  print(summary(combined_data))
  
  return(combined_data)
}

# ============================================================================
# GEMMA Testing Functions
# ============================================================================

#' Run default GEMMA on a single trait-SNP pair
#'
#' @param phenotype Vector of phenotype values
#' @param snp_genotype Matrix with single SNP genotype (n × 1)
#' @param kinship Kinship matrix
#' @param covariates Covariate matrix (default intercept only)
#' @return P-value from Wald test
run_default_gemma <- function(phenotype, snp_genotype, kinship, 
                              covariates = NULL) {
  
  # Prepare covariates (intercept if not provided)
  if (is.null(covariates)) {
    W <- matrix(1, nrow = length(phenotype), ncol = 1)
  } else {
    W <- as.matrix(covariates)
  }
  
  # Prepare data
  data <- prepare_gemma_data(
    pheno = phenotype,
    geno = snp_genotype,
    covars = W,
    kinship = kinship
  )
  
  # Run GEMMA
  result <- tryCatch({
    pygemma(
      Y = data$Y,
      X = data$X,
      W = data$W,
      K = data$K,
      verbose = 0,
      grid = FALSE,
      nproc = 1
    )
  }, error = function(e) {
    return(data.frame(p_wald = NA_real_))
  })
  
  return(result$p_wald[1])
}

#' Run modified (whitened) GEMMA on a single trait-SNP pair
#'
#' @param phenotype Vector of phenotype values
#' @param snp_genotype Matrix with single SNP genotype (n × 1)
#' @param kinship Kinship matrix
#' @param phenotype_variance Vector of phenotype variances for weighting
#' @param lambda_shrink Shrinkage parameter for variance (default 0.2)
#' @param covariates Covariate matrix (default intercept only)
#' @return P-value from Wald test
run_modified_gemma <- function(phenotype, snp_genotype, kinship,
                               phenotype_variance,
                               lambda_shrink = 0,
                               covariates = NULL) {
  
  # If no variance provided, use sample variance
  if (is.null(phenotype_variance)) {
    phenotype_variance <- rep(var(phenotype), length(phenotype))
  }
  
  # Shrink variance estimates
  log_s2 <- log(phenotype_variance)
  mu <- mean(log_s2)
  log_s2_shr <- (1 - lambda_shrink) * log_s2 + lambda_shrink * mu
  s2_shr <- exp(log_s2_shr)
  
  # Create weighting matrix
  w <- 1 / sqrt(s2_shr)
  Ww <- diag(as.vector(w))
  
  # Whiten all data
  K_whitened <- Ww %*% as.matrix(kinship) %*% Ww
  Y_w <- Ww %*% phenotype
  X_w <- Ww %*% snp_genotype
  
  # Prepare covariates
  if (is.null(covariates)) {
    W_fe <- matrix(1, nrow = length(phenotype), ncol = 1)
  } else {
    W_fe <- as.matrix(covariates)
  }
  W_fe_w <- Ww %*% W_fe
  
  # Prepare data
  data <- prepare_gemma_data(
    pheno = Y_w,
    geno = X_w,
    covars = W_fe_w,
    kinship = K_whitened
  )
  
  # Run GEMMA
  result <- tryCatch({
    pygemma(
      Y = data$Y,
      X = data$X,
      W = W_fe_w,
      K = data$K,
      verbose = 0,
      grid = FALSE,
      nproc = 1
    )
  }, error = function(e) {
    return(data.frame(p_wald = NA_real_))
  })
  
  return(result$p_wald[1])
}

#' Compile power simulation data (similar to null compilation)
#'
#' @param base_dir Base directory containing simulation files
#' @param h2_levels Vector of h2 values to process (e.g., c(0.005, 0.01, 0.05))
#' @param n_sims_per_h2 Number of simulations per h2 level
#' @param file_pattern_pheno Pattern for phenotype files (use %s for h2_int, %d for sim_id)
#' @param file_pattern_se2 Pattern for variance files
#' @param strain_names Vector of strain names
#' @param n_cores Number of cores for parallel processing
#' @return Data frame with all power simulation data
#' @examples
#' # For files like: PhenotypeSimulated_singleSNP_h20005_UHetero1sc_sampledSNP_1.txt
#' # Use pattern: "PhenotypeSimulated_singleSNP_h2%s_UHetero1sc_sampledSNP_%d.txt"
compile_power_simulations <- function(base_dir, h2_levels, n_sims_per_h2,
                                      file_pattern_pheno, file_pattern_se2,
                                      strain_names = NULL, n_cores = 1,
                                      include_null = FALSE,
                                      null_pattern_pheno = NULL,
                                      null_pattern_se2 = NULL) {
  
  cat("=== Compiling Power Simulation Data ===\n\n")
  cat(sprintf("H2 levels: %s\n", paste(h2_levels, collapse = ", ")))
  if (include_null) cat("Including null simulations (h2 = 0)\n")
  cat(sprintf("Simulations per H2: %d\n", n_sims_per_h2))
  cat(sprintf("Total files to read: %d\n", length(h2_levels) * n_sims_per_h2 * 2))
  
  all_data <- list()
  
  # ---- Null case (no h2 parameter in filename) ----
  if (include_null) {
    if (is.null(null_pattern_pheno) || is.null(null_pattern_se2)) {
      stop("Must provide null_pattern_pheno and null_pattern_se2 when include_null = TRUE")
    }
    
    cat("\nProcessing null simulations (h2 = 0)...\n")
    
    null_data <- lapply(1:n_sims_per_h2, function(sim_id) {
      pheno_file <- file.path(base_dir, sprintf(null_pattern_pheno, sim_id))
      se2_file   <- file.path(base_dir, sprintf(null_pattern_se2,   sim_id))
      
      if (!file.exists(pheno_file) || !file.exists(se2_file)) {
        if (sim_id == 1) {
          cat(sprintf("  Warning: Null files not found for sim_id=%d\n", sim_id))
          cat(sprintf("    Pheno: %s (exists: %s)\n", pheno_file, file.exists(pheno_file)))
          cat(sprintf("    SE2:   %s (exists: %s)\n", se2_file,   file.exists(se2_file)))
        }
        return(NULL)
      }
      
      tryCatch({
        pheno <- read.table(pheno_file, header = FALSE, stringsAsFactors = FALSE)
        se2   <- read.table(se2_file,   header = FALSE, stringsAsFactors = FALSE)
        
        if (ncol(pheno) >= 2 && ncol(se2) >= 2) {
          data <- data.frame(strain = pheno[, 1],
                             value  = as.numeric(pheno[, 2]),
                             se2    = as.numeric(se2[, 2]),
                             stringsAsFactors = FALSE)
        } else if (ncol(pheno) >= 2) {
          data <- data.frame(strain = pheno[, 1],
                             value  = as.numeric(pheno[, 2]),
                             se2    = as.numeric(se2[, 1]),
                             stringsAsFactors = FALSE)
        } else {
          data <- data.frame(
            strain = if (!is.null(strain_names)) strain_names else 1:nrow(pheno),
            value  = as.numeric(pheno[, 1]),
            se2    = as.numeric(se2[, 1]),
            stringsAsFactors = FALSE
          )
        }
        
        data$h2     <- 0
        data$h2_int <- "null"
        data$sim_id <- sim_id
        data$trait  <- sprintf("Trait_null_sim_%d", sim_id)
        
        return(data)
      }, error = function(e) {
        if (sim_id <= 5) cat(sprintf("  Error reading null sim_id=%d: %s\n", sim_id, e$message))
        return(NULL)
      })
    })
    
    null_data <- null_data[!sapply(null_data, is.null)]
    if (length(null_data) > 0) {
      all_data[["null"]] <- bind_rows(null_data)
      cat(sprintf("  Successfully read %d/%d null simulations\n",
                  length(null_data), n_sims_per_h2))
    } else {
      cat("  WARNING: No null files successfully read\n")
    }
  }
  
  # ---- Non-null h2 levels (existing logic unchanged) ----
  for (h2 in h2_levels) {
    h2_int <- h2
    cat(sprintf("\nProcessing h2 = %.4f (file code: h2%s)...\n", h2, h2_int))
    
    h2_data <- lapply(1:n_sims_per_h2, function(sim_id) {
      pheno_file <- file.path(base_dir, sprintf(file_pattern_pheno, h2_int, sim_id))
      se2_file   <- file.path(base_dir, sprintf(file_pattern_se2,   h2_int, sim_id))
      
      if (!file.exists(pheno_file) || !file.exists(se2_file)) {
        if (sim_id == 1) {
          cat(sprintf("  Warning: Files not found for sim_id=%d\n", sim_id))
          cat(sprintf("    Pheno: %s (exists: %s)\n", pheno_file, file.exists(pheno_file)))
          cat(sprintf("    SE2:   %s (exists: %s)\n", se2_file,   file.exists(se2_file)))
        }
        return(NULL)
      }
      
      tryCatch({
        pheno <- read.table(pheno_file, header = FALSE, stringsAsFactors = FALSE)
        se2   <- read.table(se2_file,   header = FALSE, stringsAsFactors = FALSE)
        
        if (ncol(pheno) >= 2 && ncol(se2) >= 2) {
          data <- data.frame(strain = pheno[, 1],
                             value  = as.numeric(pheno[, 2]),
                             se2    = as.numeric(se2[, 2]),
                             stringsAsFactors = FALSE)
        } else if (ncol(pheno) >= 2) {
          data <- data.frame(strain = pheno[, 1],
                             value  = as.numeric(pheno[, 2]),
                             se2    = as.numeric(se2[, 1]),
                             stringsAsFactors = FALSE)
        } else if (ncol(pheno) == 1 && ncol(se2) == 1) {
          data <- data.frame(
            strain = if (!is.null(strain_names)) strain_names else 1:nrow(pheno),
            value  = as.numeric(pheno[, 1]),
            se2    = as.numeric(se2[, 1]),
            stringsAsFactors = FALSE
          )
        } else {
          data <- data.frame(
            strain = if (!is.null(strain_names)) strain_names else 1:nrow(pheno),
            value  = as.numeric(pheno[, ncol(pheno)]),
            se2    = as.numeric(se2[, ncol(se2)]),
            stringsAsFactors = FALSE
          )
        }
        
        data$h2     <- h2
        data$h2_int <- h2_int
        data$sim_id <- sim_id
        data$trait  <- sprintf("Trait_h2_%s_sim_%d", h2_int, sim_id)
        
        return(data)
      }, error = function(e) {
        if (sim_id <= 5) cat(sprintf("  Error reading sim_id=%d: %s\n", sim_id, e$message))
        return(NULL)
      })
    })
    
    h2_data <- h2_data[!sapply(h2_data, is.null)]
    if (length(h2_data) > 0) {
      all_data[[as.character(h2)]] <- bind_rows(h2_data)
      cat(sprintf("  Successfully read %d/%d simulations\n",
                  length(h2_data), n_sims_per_h2))
    } else {
      cat(sprintf("  WARNING: No files successfully read for h2 = %.4f\n", h2))
    }
  }
  
  # ---- Combine and summarise ----
  combined_data <- bind_rows(all_data)
  
  if (nrow(combined_data) == 0) {
    stop("No data successfully compiled! Check file paths and patterns.")
  }
  
  cat(sprintf("\nTotal observations compiled: %s\n",
              format(nrow(combined_data), big.mark = ",")))
  cat(sprintf("Unique h2 levels: %s\n",
              paste(sort(unique(combined_data$h2)), collapse = ", ")))
  cat(sprintf("Unique strains: %d\n", length(unique(combined_data$strain))))
  
  return(combined_data)
}

# ============================================================================
# Power Calculation Functions
# ============================================================================

#' Calculate power for both GEMMA methods
#'
#' @param power_data Data frame with power simulations
#' @param genotype_matrix Full genotype matrix
#' @param kinship Kinship matrix
#' @param causal_snp_map Data frame mapping sim_id to causal_snp
#' @param alpha Significance threshold
#' @param n_cores Number of cores for parallel processing
#' @return Data frame with power test results
calculate_power_both_methods <- function(power_data, genotype_matrix, kinship,
                                         causal_snp_map = NULL, alpha = 0.05,
                                         n_cores = 1) {
  
  cat("=== Calculating Power for Both Methods ===\n\n")
  
  # Get unique trait-h2-sim combinations
  test_combos <- power_data %>%
    distinct(trait, h2, sim_id) %>%
    arrange(h2, sim_id)
  
  n_tests <- nrow(test_combos)
  cat(sprintf("Running %d power tests...\n", n_tests))
  
  # Source GEMMA functions (assuming they're available)
  # Ensure run_default_gemma and run_modified_gemma are loaded
  
  # Function to test a single trait
  test_single_trait <- function(i) {
    h2 <- test_combos$h2[i]
    sim_id <- test_combos$sim_id[i]
    trait_name <- test_combos$trait[i]
    
    # Get phenotype data
    trait_data <- power_data %>%
      filter(trait == trait_name) %>%
      arrange(strain)
    
    # Match strain order to genotype matrix
    if (!all(trait_data$strain %in% rownames(genotype_matrix))) {
      return(NULL)
    }
    
    trait_data <- trait_data %>%
      arrange(match(strain, rownames(genotype_matrix)))
    
    phenotype <- trait_data$value
    phenotype_var <- trait_data$se2
    
    # Get causal SNP
    if (!is.null(causal_snp_map)) {
      causal_snp <- causal_snp_map %>%
        filter(sim_id == !!sim_id) %>%
        pull(causal_snp)
      
      if (length(causal_snp) == 0) {
        return(NULL)
      }
    } else {
      # Assume column name pattern or use first SNP
      causal_snp <- colnames(genotype_matrix)[sim_id]
    }
    
    if (!causal_snp %in% colnames(genotype_matrix)) {
      return(NULL)
    }
    
    snp_geno <- genotype_matrix[, causal_snp, drop = FALSE]
    
    # Run both methods
    p_default <- run_default_gemma(phenotype, snp_geno, kinship)
    p_modified <- run_modified_gemma(phenotype, snp_geno, kinship,
                                     phenotype_variance = phenotype_var)
    
    data.frame(
      h2 = h2,
      sim_id = sim_id,
      trait = trait_name,
      causal_snp = causal_snp,
      p_default = p_default,
      p_modified = p_modified,
      detected_default = p_default < alpha,
      detected_modified = p_modified < alpha,
      log10p_default = -log10(p_default),
      log10p_modified = -log10(p_modified),
      log10p_ratio = -log10(p_modified) / -log10(p_default),
      stringsAsFactors = FALSE
    )
  }
  
  # Run tests (parallel or sequential)
  if (n_cores > 1) {
    cat(sprintf("Using %d cores...\n", n_cores))
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)
    
    clusterExport(cl, c("power_data", "genotype_matrix", "kinship",
                        "causal_snp_map", "alpha", "test_combos",
                        "run_default_gemma", "run_modified_gemma",
                        "prepare_gemma_data", "pygemma",
                        "calc_lambda_restricted", "calc_beta_vg_ve_restricted",
                        "likelihood_derivative1_lambda", 
                        "likelihood_derivative2_lambda",
                        "likelihood_lambda"),
                  envir = environment())
    
    clusterEvalQ(cl, library(tidyverse))
    
    results_list <- pbapply::pblapply(1:n_tests, test_single_trait, cl = cl)
  } else {
    if (requireNamespace("pbapply", quietly = TRUE)) {
      results_list <- pbapply::pblapply(1:n_tests, test_single_trait)
    } else {
      results_list <- lapply(1:n_tests, function(i) {
        if (i %% 50 == 0) cat(sprintf("  Progress: %d/%d\n", i, n_tests))
        test_single_trait(i)
      })
    }
  }
  
  # Combine results
  results_df <- bind_rows(results_list)
  
  cat("\nPower calculation complete!\n")
  
  return(results_df)
}

# Heritability levels tested (as decimals)
#h2_levels <- c("0005", "001", "0015", "002", "0025", "003", "0035", "004", "0045", "005")  # Adjust to your actual levels
h2_levels <- NA  # Adjust to your actual levels

# Number of simulations per h2 level
n_sims_per_h2 <- 1000

# Base directory containing simulation files
base_dir <- "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping"

# File patterns using %s for h2_int (e.g., "0005") and %d for sim_id
# Example file: PhenotypeSimulated_singleSNP_h20005_UHetero1sc_sampledSNP_1.txt

# Output directory
output_dir <- "/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/QTL_mapping/fdr_analysis"

# Significance threshold
alpha <- 0.05

# Number of cores
n_cores <- 8

# ============================================================================
# Load Required Data
# ============================================================================

cat("=== Loading Required Data ===\n\n")

# Load genotype matrix
cat("Loading genotype matrix...\n")
load(file = "NORMgenomat.samp.Rda")  # Contains 'g'
g <- NORMAL.genomat.samp_mat
cat(sprintf("  Loaded: %d strains × %d SNPs\n", nrow(g), ncol(g)))

# Load kinship matrix
cat("Loading kinship matrix...\n")
fread(file = "gemmaGRM.full.sXX.txt") -> kinship
kinship_matrix <- as.matrix(kinship)
# kinship_matrix_PD <- nearPD(kinship_matrix, corr = FALSE,
#                           keepDiag = TRUE,
#                           ensureSymmetry = TRUE)$mat
cat(sprintf("  Loaded: %d × %d\n", nrow(kinship_matrix), ncol(kinship_matrix)))

load("rlnorm_mneg3sd4_simulations_compiled_060226.Rda")
load("LogNorm_mneg15sd25_simulations_compiled_052826.Rda")
load("./fdr_analysis/FITmeanvar_null_simulations_compiled_052826.Rda")
load("./fdr_analysis/UHetero1sc_null_simulations_compiled_052826.Rda")

plotting_data <- power_data %>% 
  filter(trait == "Trait_null_sim_1")
plot(x = plotting_data$value, y = plotting_data$se2, 
     xlab = "Average Simulated Trait Value per Strain", ylab = "Fitted Variance in Simulated Trait Value per Strain",
     main = "Mean-Variance Relationship for Simulated Null Trait Data, sim 1"
     )
hist(x = plotting_data$value,
     xlab = "Average Simulated Trait Value per Strain",
     main = "Simulated Null Trait Values, sim 1")
hist(x = plotting_data$se2,
     xlab = "Variance in Simulated Trait Value per Strain",
     main = "Simulated Null Variances drawn from Log-Normal distribution, sim 1")

loess_fit <- loess(se2 ~ value, data = plotting_data, span = 0.80)
x_seq     <- seq(min(plotting_data$value), max(plotting_data$value), length.out = 200)

pred      <- predict(loess_fit, newdata = data.frame(value = x_seq), se = TRUE)
ci_upper  <- pred$fit + 1.96 * pred$se.fit
ci_lower  <- pred$fit - 1.96 * pred$se.fit

plot(x = plotting_data$value, y = plotting_data$se2,
     xlab = "Average Simulated Trait Value per Strain",
     ylab = "Fitted Variance in Simulated Trait Value per Strain",
     main = "Mean-Variance Relationship for Simulated Null Trait Data, sim 1"
)

lines(x_seq, pred$fit, col = "firebrick", lwd = 2)
lines(x_seq, ci_upper, col = "firebrick", lwd = 1, lty = 2)
lines(x_seq, ci_lower, col = "firebrick", lwd = 1, lty = 2)

# ggplot2 ver
plotting_data <- power_data %>%
  filter(trait == "Trait_null_sim_1")

# 1. Scatter plot: Mean-Variance Relationship
p1 <- ggplot(plotting_data, aes(x = value, y = se2)) +
  geom_point(size = 2) +
  labs(
    x     = "Average Simulated Trait Value per Strain",
    #y     = "Fitted Variance in Simulated Trait Value per Strain",
    y     = "Variance in Simulated Trait Value per Strain",
    title = "Mean-Variance Relationship for Simulated Null Trait Data with Lognorm(-3,4) variances, sim 1"
  )   +
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
p1

# 2. Histogram of trait values
p2 <- ggplot(plotting_data, aes(x = value)) +
  geom_histogram(bins = 30) +
  labs(
    x     = "Average Simulated Trait Value per Strain",
    title = "Simulated Null Trait Values, sim 1"
  ) +
  theme_bw()

# 3. Histogram of variances
p3 <- ggplot(plotting_data, aes(x = se2)) +
  geom_histogram(bins = 30) +
  labs(
    x     = "Variance in Simulated Trait Value per Strain",
    title = "Simulated Null Variances drawn from Log-Normal distribution, sim 1"
  ) +
  theme_bw()

# 4. Scatter plot with LOESS smoother + 95% CI
p4 <- ggplot(plotting_data, aes(x = value, y = se2)) +
  geom_point() +
  geom_smooth(
    method  = "loess",
    span    = 0.80,
    color   = "firebrick",
    fill    = "firebrick",
    alpha   = 0.15,         # shaded CI ribbon
    linewidth = 1
  ) +
  labs(
    x     = "Average Simulated Trait Value per Strain",
    y     = "Fitted Variance in Simulated Trait Value per Strain",
    title = "Mean-Variance Relationship for Simulated Null Trait Data, sim 1"
  ) +
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

# Print all four plots
print(p1)
print(p2)
print(p3)
print(p4)

# ============================================================================
# Run Power Tests for Both Methods
# ============================================================================

cat("\n=== Running Power Tests ===\n\n")

# Check if results already exist
results_file <- file.path(output_dir, "LogNorm_mneg15sd25_nullsims_results_0602.Rda")
results_file <- file.path(output_dir, "FITmeanvar_null_results_060126.Rda")

# Null Homogenous
power_data$se2 <- 1
results_file <- file.path(output_dir, "UHetero_null_homog_results_060126.Rda")

causal_snp_map <- NULL

if (file.exists(results_file)) {
  cat("Loading existing power results...\n")
  load(results_file)
  cat(sprintf("Loaded %d test results\n", nrow(power_results)))
  
} else {
  cat("Calculating power for both methods...\n")
  
  power_results <- calculate_power_both_methods(
    power_data = power_data,
    genotype_matrix = g,
    kinship = kinship_matrix,
    #kinship = diag(x = 1, nrow = 192, ncol = 192),
    causal_snp_map = causal_snp_map,
    alpha = alpha,
    n_cores = n_cores
  )
  
  # Save results
  save(power_results, file = results_file)
  cat(sprintf("\nSaved results to: %s\n", results_file))
}

load("./fdr_analysis/LogNorm_mneg15sd25_nullsims_results_0602.Rda")
load("./fdr_analysis/LogNorm_mneg3sd4_nullsims_results_0602.Rda")
load("./fdr_analysis/UHetero1sc_null_results_052826.Rda")
load("./fdr_analysis/FITmeanvar_null_results_052826.Rda")

power_results %>% 
  summarise(sum(p_default < 0.05) / 1000,
            sum(p_modified < 0.05) / 1000)

mc_table <- table(
  default  = power_results$detected_default,
  modified = power_results$detected_modified
)
print(mc_table)

mcnemar.test(mc_table)

# Filter to null simulations
null_data <- power_results

median(qchisq(1 - null_data$p_default,  1)) 
median(null_data$p_default)

median(qchisq(1 - null_data$p_modified,  1)) 
median(null_data$p_modified)

# Compute genomic inflation factor λ
lambda_df <- tibble(
  model  = c("GEMMA", "GEMMAmod"),
  lambda = c(
    median(qchisq(1 - null_data$p_default,  1), na.rm = TRUE) / 0.4549,
    median(qchisq(1 - null_data$p_modified, 1), na.rm = TRUE) / 0.4549
  )
) %>%
  mutate(label = paste0("λ = ", round(lambda, 3)))
lambda_df

library(gap)
lambda_default  <- gc.lambda(null_data$p_default)
lambda_modified <- gc.lambda(null_data$p_modified)
lambda_default
lambda_modified

hist(null_data$p_default,
     main = "Null p-values for GEMMA: Var~LogNorm(-1.5,2.5)",
     xlab = "p-values",
     ylim = c(0,200))
hist(null_data$p_modified,
     main = "Null p-values for GEMMAvar: Var~LogNorm(-1.5,2.5)",
     xlab = "p-values",
     ylim = c(0,200))

# Build long-format expected vs observed
qq_df <- bind_rows(
  null_data %>%
    arrange(p_default) %>%
    mutate(
      model = "GEMMA",
      expected = -log10(ppoints(n())),
      observed = -log10(p_default)
    ),
  null_data %>%
    arrange(p_modified) %>%
    mutate(
      model = "GEMMAvar",
      expected = -log10(ppoints(n())),
      observed = -log10(p_modified)
    )
) %>%
  select(model, expected, observed)

# Determine label position (upper-left of each facet)
label_x <- min(qq_df$expected) + 0.05 * diff(range(qq_df$expected))
label_y_GEMMA <- max(qq_df$observed) * 0.95
label_y_GEMMAmod <- max(qq_df$observed) * 0.9

# Plot
ggplot(qq_df, aes(x = expected, y = observed, color = model)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +
  # geom_text(
  #   data = lambda_df[lambda_df$model == "GEMMA",],
  #   aes(x = label_x, y = label_y_GEMMA, label = label),
  #   inherit.aes = FALSE,
  #   hjust = 0, vjust = 1,
  #   size = 4, fontface = "bold",
  #   color = c("#d6604d")) + 
  # geom_text(
  #   data = lambda_df[lambda_df$model == "GEMMAmod",],
  #   aes(x = label_x, y = label_y_GEMMAmod, label = label),
  #   inherit.aes = FALSE,
  #   hjust = 0, vjust = 1,
  #   size = 4, fontface = "bold",
  #   color = c("#2166ac")) + 
  scale_color_manual(
    values = c("GEMMA" = "#d6604d", "GEMMAvar" = "#2166ac" )
  ) +
  labs(
    #title    = "QQ-plot: GEMMA vs GEMMAvar (null simulations, Var~LogNormal(-1.5,2.5))",
    title    = "QQ-plot: GEMMA vs GEMMAvar (null simulations, Var~LogNormal(-3,4))",
    #title    = "QQ-plot: GEMMA vs GEMMAvar (null simulations, Var~U(0,1.5))",
    #title    = "QQ-plot: GEMMA vs GEMMAvar (null simulations, Fitted mean-variance)",
    x        = expression(Expected ~ -log[10](p)),
    y        = expression(Observed ~ -log[10](p)),
    color    = "Model"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggplot(data = power_results, aes(x = log10p_default, y = log10p_modified)) + geom_point() +
  #labs(y = "GEMMAvar -log10(p-value)", x = "GEMMA -log10(p-value)", title = "Simulated Null Normally Distributed traits with homogenous variances and 1000 different SNPs") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated Negative Binomial with random strain effects, dispersion=0.04-0.05, logFC= 0.3 mean-var noise traits with 1000 different SNPs") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated Gamma(B = 0.25, s = 4) traits, heterogenous variance p-values using Garter Snake genotypes") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated h2 = 0.05 traits, imposed mean-var relationship p-values using Garter Snake genotypes") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated h2 = 0.025 traits, fitted mean-var relationship p-values using Garter Snake genotypes") +
  #labs(y = "Modified R ver. -log10p", x = "Default R ver. -log10p", title = "Simulated h2 = 0.015 traits, fitted mean-var relationship p-values using Garter Snake genotypes") +
  geom_abline(intercept=0, slope=1, color = "red") + 
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

