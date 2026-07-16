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

# # Compile all data
simulated_null_data <- compile_simulated_data(
  n_sims = 1000,
  pheno_template = pheno_template,
  se2_template = se2_template,
  strain_names = NULL,  # Will try to extract from files
  n_cores = 8,
  sample_first_file = TRUE
)

metadata <- data.frame(
  sim_id = c(1:1000),
  causal_snp = colnames(NORMAL.genomat.samp_mat)
)

simulated_null_data <- simulated_null_data %>% 
  select(!causal_snp) %>% 
  left_join(metadata)

strain_col <- rep(rownames(NORMAL.genomat.samp_mat), 1000)

simulated_null_data$strain <- strain_col

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


# ============================================================================
# Power Analysis for GEMMA Methods
# ============================================================================
# Functions to calculate and compare power between Default and Modified GEMMA
# across different heritability levels
# ============================================================================

library(tidyverse)
library(parallel)

# ============================================================================
# Data Compilation for Power Simulations
# ============================================================================

#' Convert h2 decimal to integer format used in filenames
#'
#' @param h2 Heritability value (e.g., 0.005, 0.01, 0.5)
#' @return Integer representation (e.g., 0005, 0010, 0500)
h2_to_int <- function(h2) {
  # Convert h2 to 4-digit integer format
  # 0.005 -> 0005, 0.01 -> 0010, 0.5 -> 0500, etc.
  h2_int <- round(h2 * 10000)
  sprintf("%04d", h2_int)
}

#' Convert integer h2 format back to decimal
#'
#' @param h2_int Integer h2 string (e.g., "0005")
#' @return Decimal h2 value (e.g., 0.0005)
int_to_h2 <- function(h2_int) {
  as.numeric(h2_int) / 10000
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
                                      strain_names = NULL, n_cores = 1) {
  
  cat("=== Compiling Power Simulation Data ===\n\n")
  cat(sprintf("H2 levels: %s\n", paste(h2_levels, collapse = ", ")))
  cat(sprintf("Simulations per H2: %d\n", n_sims_per_h2))
  cat(sprintf("Total files to read: %d\n", length(h2_levels) * n_sims_per_h2 * 2))
  
  all_data <- list()
  
  for (h2 in h2_levels) {
    #cat(sprintf("\nProcessing h2 = %.4f...\n", h2))
    
    # Convert h2 to integer format for filename
    h2_int <- h2
    cat(sprintf("  File code: h2%s\n", h2_int))
    
    # Read simulations for this h2 level
    h2_data <- lapply(1:n_sims_per_h2, function(sim_id) {
      # Build file paths with h2_int as string and sim_id as integer
      pheno_file <- file.path(base_dir, sprintf(file_pattern_pheno, h2_int, sim_id))
      se2_file <- file.path(base_dir, sprintf(file_pattern_se2, h2_int, sim_id))
      
      if (!file.exists(pheno_file) || !file.exists(se2_file)) {
        if (sim_id == 1) {
          # Report first missing file for debugging
          cat(sprintf("  Warning: Files not found for sim_id=%d\n", sim_id))
          cat(sprintf("    Pheno: %s (exists: %s)\n", pheno_file, file.exists(pheno_file)))
          cat(sprintf("    SE2:   %s (exists: %s)\n", se2_file, file.exists(se2_file)))
        }
        return(NULL)
      }
      
      tryCatch({
        pheno <- read.table(pheno_file, header = FALSE, stringsAsFactors = FALSE)
        se2 <- read.table(se2_file, header = FALSE, stringsAsFactors = FALSE)
        
        # Combine - handle different file formats
        if (ncol(pheno) >= 2 && ncol(se2) >= 2) {
          # Files have strain IDs in first column
          data <- data.frame(
            strain = pheno[, 1],
            value = as.numeric(pheno[, 2]),
            se2 = as.numeric(se2[, 2]),
            stringsAsFactors = FALSE
          )
        } else if (ncol(pheno) >= 2) {
          # Pheno has strain IDs, se2 doesn't
          data <- data.frame(
            strain = pheno[, 1],
            value = as.numeric(pheno[, 2]),
            se2 = as.numeric(se2[, 1]),
            stringsAsFactors = FALSE
          )
        } else if (ncol(pheno) == 1 && ncol(se2) == 1) {
          # Neither has strain IDs
          data <- data.frame(
            strain = if (!is.null(strain_names)) strain_names else 1:nrow(pheno),
            value = as.numeric(pheno[, 1]),
            se2 = as.numeric(se2[, 1]),
            stringsAsFactors = FALSE
          )
        } else {
          # Use last column for multi-column files
          data <- data.frame(
            strain = if (!is.null(strain_names)) strain_names else 1:nrow(pheno),
            value = as.numeric(pheno[, ncol(pheno)]),
            se2 = as.numeric(se2[, ncol(se2)]),
            stringsAsFactors = FALSE
          )
        }
        
        data$h2 <- h2
        data$h2_int <- h2_int
        data$sim_id <- sim_id
        data$trait <- sprintf("Trait_h2_%s_sim_%d", h2_int, sim_id)
        
        return(data)
      }, error = function(e) {
        if (sim_id <= 5) {
          cat(sprintf("  Error reading sim_id=%d: %s\n", sim_id, e$message))
        }
        return(NULL)
      })
    })
    
    # Remove NULLs and combine
    h2_data <- h2_data[!sapply(h2_data, is.null)]
    if (length(h2_data) > 0) {
      all_data[[as.character(h2)]] <- bind_rows(h2_data)
      cat(sprintf("  Successfully read %d/%d simulations\n", 
                  length(h2_data), n_sims_per_h2))
    } else {
      cat(sprintf("  WARNING: No files successfully read for h2 = %.4f\n", h2))
    }
  }
  
  # Combine all h2 levels
  combined_data <- bind_rows(all_data)
  
  if (nrow(combined_data) == 0) {
    stop("No data successfully compiled! Check file paths and patterns.")
  }
  
  cat(sprintf("\nTotal observations compiled: %s\n", 
              format(nrow(combined_data), big.mark = ",")))
  cat(sprintf("Unique h2 levels: %d\n", length(unique(combined_data$h2))))
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

#' Summarize power by h2 level
#'
#' @param power_results Results from calculate_power_both_methods
#' @param alpha Significance threshold
#' @return Data frame with power summary by h2
summarize_power_by_h2 <- function(power_results, alpha = 0.05) {
  
  power_summary <- power_results %>%
    group_by(h2) %>%
    summarise(
      n_tests = n(),
      
      # Power (detection rate)
      power_default = mean(detected_default, na.rm = TRUE),
      power_modified = mean(detected_modified, na.rm = TRUE),
      power_ratio = power_modified / power_default,
      power_difference = power_modified - power_default,
      
      # Both methods detect
      both_detected = mean(detected_default & detected_modified, na.rm = TRUE),
      
      # Only one method detects
      only_default = mean(detected_default & !detected_modified, na.rm = TRUE),
      only_modified = mean(!detected_default & detected_modified, na.rm = TRUE),
      
      # Neither detects
      neither_detected = mean(!detected_default & !detected_modified, na.rm = TRUE),
      
      # P-value comparisons (among detected signals)
      mean_log10p_default = mean(log10p_default[detected_default], na.rm = TRUE),
      mean_log10p_modified = mean(log10p_modified[detected_modified], na.rm = TRUE),
      
      # P-value ratio (how much stronger modified is)
      mean_log10p_ratio = mean(log10p_ratio[detected_default & detected_modified], 
                               na.rm = TRUE),
      median_log10p_ratio = median(log10p_ratio[detected_default & detected_modified], 
                                   na.rm = TRUE),
      
      # Standard errors
      se_power_default = sqrt(power_default * (1 - power_default) / n_tests),
      se_power_modified = sqrt(power_modified * (1 - power_modified) / n_tests),
      
      .groups = 'drop'
    ) %>%
    mutate(
      power_ratio_lower_ci = power_ratio - 1.96 * se_power_modified / power_default,
      power_ratio_upper_ci = power_ratio + 1.96 * se_power_modified / power_default
    )
  
  return(power_summary)
}

#' Calculate power gain statistics
#'
#' @param power_results Results from calculate_power_both_methods
#' @return List with detailed power gain statistics
calculate_power_gain_statistics <- function(power_results) {
  
  cat("\n=== Power Gain Statistics ===\n\n")
  
  # Overall statistics across all h2 levels
  overall_stats <- power_results %>%
    summarise(
      total_tests = n(),
      
      # Power across all h2
      overall_power_default = mean(detected_default, na.rm = TRUE),
      overall_power_modified = mean(detected_modified, na.rm = TRUE),
      
      # Concordance
      both_detected = sum(detected_default & detected_modified, na.rm = TRUE),
      only_default = sum(detected_default & !detected_modified, na.rm = TRUE),
      only_modified = sum(!detected_default & detected_modified, na.rm = TRUE),
      neither = sum(!detected_default & !detected_modified, na.rm = TRUE),
      
      # Percentage gains
      pct_both = 100 * both_detected / total_tests,
      pct_only_default = 100 * only_default / total_tests,
      pct_only_modified = 100 * only_modified / total_tests,
      pct_neither = 100 * neither / total_tests
    )
  
  cat("Overall Power:\n")
  cat(sprintf("  Default GEMMA:  %.2f%%\n", 
              100 * overall_stats$overall_power_default))
  cat(sprintf("  Modified GEMMA: %.2f%%\n", 
              100 * overall_stats$overall_power_modified))
  cat(sprintf("  Absolute gain:  %.2f%%\n",
              100 * (overall_stats$overall_power_modified - 
                       overall_stats$overall_power_default)))
  cat(sprintf("  Relative gain:  %.2f%%\n\n",
              100 * (overall_stats$overall_power_modified / 
                       overall_stats$overall_power_default - 1)))
  
  cat("Detection Patterns:\n")
  cat(sprintf("  Both methods:     %5d (%.1f%%)\n", 
              overall_stats$both_detected, overall_stats$pct_both))
  cat(sprintf("  Default only:     %5d (%.1f%%)\n", 
              overall_stats$only_default, overall_stats$pct_only_default))
  cat(sprintf("  Modified only:    %5d (%.1f%%)\n", 
              overall_stats$only_modified, overall_stats$pct_only_modified))
  cat(sprintf("  Neither:          %5d (%.1f%%)\n\n", 
              overall_stats$neither, overall_stats$pct_neither))
  
  # P-value significance gain
  pvalue_stats <- power_results %>%
    filter(detected_default & detected_modified) %>%
    summarise(
      n_both_detected = n(),
      
      # How much more significant is modified?
      mean_p_fold_change = mean(p_default / p_modified, na.rm = TRUE),
      median_p_fold_change = median(p_default / p_modified, na.rm = TRUE),
      
      # -log10(p) differences
      mean_log10p_gain = mean(log10p_modified - log10p_default, na.rm = TRUE),
      median_log10p_gain = median(log10p_modified - log10p_default, na.rm = TRUE),
      
      # Percentage with modified > default
      pct_modified_stronger = 100 * mean(log10p_modified > log10p_default, na.rm = TRUE)
    )
  
  cat("P-value Significance (among signals detected by both):\n")
  cat(sprintf("  N = %d\n", pvalue_stats$n_both_detected))
  cat(sprintf("  Mean p-value fold change:   %.2fx smaller (modified vs default)\n", 
              pvalue_stats$mean_p_fold_change))
  cat(sprintf("  Median p-value fold change: %.2fx smaller\n", 
              pvalue_stats$median_p_fold_change))
  cat(sprintf("  Mean -log10(p) gain:        %.2f\n", 
              pvalue_stats$mean_log10p_gain))
  cat(sprintf("  Median -log10(p) gain:      %.2f\n", 
              pvalue_stats$median_log10p_gain))
  cat(sprintf("  Modified stronger:          %.1f%% of cases\n\n", 
              pvalue_stats$pct_modified_stronger))
  
  # By h2 level
  h2_stats <- power_results %>%
    group_by(h2) %>%
    summarise(
      n = n(),
      power_default = mean(detected_default, na.rm = TRUE),
      power_modified = mean(detected_modified, na.rm = TRUE),
      power_gain = power_modified - power_default,
      
      # Among both detected
      n_both = sum(detected_default & detected_modified, na.rm = TRUE),
      mean_log10p_gain = mean(log10p_modified[detected_default & detected_modified] - 
                                log10p_default[detected_default & detected_modified], 
                              na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("Power by H2 Level:\n")
  print(h2_stats)
  
  return(list(
    overall = overall_stats,
    pvalue_stats = pvalue_stats,
    by_h2 = h2_stats
  ))
}

# ============================================================================
# Visualization Functions
# ============================================================================

#' Plot power curves for both methods
#'
#' @param power_summary Summary from summarize_power_by_h2
#' @return ggplot object
plot_power_curves <- function(power_summary) {
  
  # Reshape for plotting
  plot_data <- power_summary %>%
    select(h2, power_default, power_modified) %>%
    pivot_longer(cols = c(power_default, power_modified),
                 names_to = "Method",
                 values_to = "Power") %>%
    mutate(Method = recode(Method,
                           power_default = "Default GEMMA",
                           power_modified = "Modified GEMMA"))
  
  p <- ggplot(plot_data, aes(x = h2, y = Power, color = Method, 
                             shape = Method, linetype = Method)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(breaks = unique(plot_data$h2)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", 
               color = "gray50", alpha = 0.5) +
    annotate("text", x = min(plot_data$h2), y = 0.8, 
             label = "80% Power", vjust = -0.5, color = "gray50", size = 3.5) +
    labs(title = "Power Curves: Default vs Modified GEMMA",
         subtitle = "Detection rate at α = 0.05",
         x = expression(Heritability~(h^2)),
         y = "Power (Detection Rate)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.position = "bottom")
  
  return(p)
}

#' Plot power ratio (Modified / Default)
#'
#' @param power_summary Summary from summarize_power_by_h2
#' @return ggplot object
plot_power_ratio <- function(power_summary) {
  
  p <- ggplot(power_summary, aes(x = h2, y = power_ratio)) +
    geom_line(size = 1.2, color = "steelblue") +
    geom_point(size = 3, color = "steelblue") +
    geom_ribbon(aes(ymin = power_ratio_lower_ci, ymax = power_ratio_upper_ci),
                alpha = 0.2, fill = "steelblue") +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = mean(power_summary$h2), y = 1.0,
             label = "Equal Power", vjust = -0.5, color = "red", size = 6) +
    scale_x_continuous(breaks = unique(power_summary$h2)) +
    labs(title = "Power Ratio: Modified / Default GEMMA",
         subtitle = "Uniform Variance",
         #subtitle = "Fitted Mean-Variance",
         x = expression(Heritability~(h^2)),
         y = "Power Ratio (Modified / Default)") +
    ylim(c(0.9, 2.4)) + 
    theme_bw()   + 
    theme(
      legend.position  = "bottom",
      axis.text.x      = element_text(size = 16),
      axis.text.y      = element_text(size = 16),
      plot.title       = element_text(hjust = 0.5, size = 22),
      plot.subtitle    = element_text(hjust = 0.5, size = 14),
      legend.title     = element_text(size = 14),
      legend.text      = element_text(size = 14),
      axis.title.x     = element_text(size = 18),
      axis.title.y     = element_text(size = 18)
    )
  
  return(p)
}

#' Plot power gain (difference)
#'
#' @param power_summary Summary from summarize_power_by_h2
#' @return ggplot object
plot_power_gain <- function(power_summary) {
  
  p <- ggplot(power_summary, aes(x = h2, y = power_difference)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks = unique(power_summary$h2)) +
    labs(title = "Power Gain: Modified - Default GEMMA",
         subtitle = "Absolute difference in detection rates",
         x = expression(Heritability~(h^2)),
         y = "Power Gain (Modified - Default)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  
  return(p)
}

#' Plot p-value significance gain
#'
#' @param power_results Full results from calculate_power_both_methods
#' @return ggplot object
plot_pvalue_gain <- function(power_results) {
  
  # Calculate mean -log10(p) gain by h2
  pvalue_gain <- power_results %>%
    filter(detected_default & detected_modified) %>%
    group_by(h2) %>%
    summarise(
      n = n(),
      mean_gain = mean(log10p_modified - log10p_default, na.rm = TRUE),
      se_gain = sd(log10p_modified - log10p_default, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  p <- ggplot(pvalue_gain, aes(x = h2, y = mean_gain)) +
    geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_gain - 1.96*se_gain, 
                      ymax = mean_gain + 1.96*se_gain),
                  width = 0.0025, size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = unique(pvalue_gain$h2)) +
    labs(title = expression(paste("P-value Significance Gain: ", 
                                  Delta, "-log"[10], "(p)")),
         subtitle = "Among signals detected by both methods (positive = Modified stronger)",
         x = expression(Heritability~(h^2)),
         y = expression(paste("Mean ", Delta, "-log"[10], "(p) = Modified - Default"))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  
  return(p)
}

#' Plot detection pattern breakdown
#'
#' @param power_summary Summary from summarize_power_by_h2
#' @return ggplot object
plot_detection_patterns <- function(power_summary) {
  
  # Reshape data
  plot_data <- power_summary %>%
    select(h2, both_detected, only_default, only_modified, neither_detected) %>%
    pivot_longer(cols = -h2, names_to = "Pattern", values_to = "Proportion") %>%
    mutate(Pattern = recode(Pattern,
                            both_detected = "Both Methods",
                            only_default = "Default Only",
                            only_modified = "Modified Only",
                            neither_detected = "Neither Method"),
           Pattern = factor(Pattern, levels = c("Both Methods", "Modified Only",
                                                "Default Only", "Neither Method")))
  
  p <- ggplot(plot_data, aes(x = h2, y = Proportion, fill = Pattern)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Both Methods" = "#2E7D32",
                                 "Modified Only" = "#1976D2",
                                 "Default Only" = "#F57C00",
                                 "Neither Method" = "#757575")) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks = unique(plot_data$h2)) +
    labs(title = "Detection Patterns by Heritability",
         x = expression(Heritability~(h^2)),
         y = "Proportion of Tests",
         fill = "Detection Pattern") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11))
  
  return(p)
}

#' Create comprehensive power analysis report
#'
#' @param power_results Full results from calculate_power_both_methods
#' @param power_summary Summary from summarize_power_by_h2
#' @param output_dir Directory to save plots
#' @return List of plots
create_power_report <- function(power_results, power_summary, output_dir) {
  
  cat("\n=== Creating Power Analysis Report ===\n\n")
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate all plots
  cat("Generating plots...\n")
  
  p1 <- plot_power_curves(power_summary)
  # ggsave(file.path(output_dir, "power_curves.png"), p1, 
  #        width = 10, height = 6, dpi = 300)
  # cat("  Saved: power_curves.png\n")
  
  p2 <- plot_power_ratio(power_summary)
  # ggsave(file.path(output_dir, "power_ratio.png"), p2, 
  #        width = 10, height = 6, dpi = 300)
  # cat("  Saved: power_ratio.png\n")
  
  p3 <- plot_power_gain(power_summary)
  # ggsave(file.path(output_dir, "power_gain.png"), p3, 
  #        width = 10, height = 6, dpi = 300)
  # cat("  Saved: power_gain.png\n")
  
  p4 <- plot_pvalue_gain(power_results)
  # ggsave(file.path(output_dir, "pvalue_significance_gain.png"), p4, 
  #        width = 10, height = 6, dpi = 300)
  # cat("  Saved: pvalue_significance_gain.png\n")
  
  p5 <- plot_detection_patterns(power_summary)
  # ggsave(file.path(output_dir, "detection_patterns.png"), p5, 
  #        width = 10, height = 6, dpi = 300)
  # cat("  Saved: detection_patterns.png\n")
  
  # Save summary tables
  # write_csv(power_summary, file.path(output_dir, "power_summary_by_h2.csv"))
  # cat("  Saved: power_summary_by_h2.csv\n")
  # 
  # write_csv(power_results, file.path(output_dir, "power_results_full.csv"))
  # cat("  Saved: power_results_full.csv\n")
  
  cat("\nReport complete!\n")
  
  return(list(
    power_curves = p1,
    power_ratio = p2,
    power_gain = p3,
    pvalue_gain = p4,
    detection_patterns = p5
  ))
}

cat("Power analysis functions loaded successfully.\n")

# ============================================================================
# Run Power Analysis - Compare Default vs Modified GEMMA
# ============================================================================

# ============================================================================
# Configuration
# ============================================================================

# Heritability levels tested (as decimals)
h2_levels <- c("0005", "001", "0015", "002", "0025", "003", "0035", "004", "0045", "005")  # Adjust to your actual levels

# Number of simulations per h2 level
n_sims_per_h2 <- 1000

# Base directory containing simulation files
base_dir <- "/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/preWhiten"

# File patterns using %s for h2_int (e.g., "0005") and %d for sim_id
# Example file: PhenotypeSimulated_singleSNP_h20005_UHetero1sc_sampledSNP_1.txt
pheno_template <- "PhenotypeSimulated_singleSNP_h2%s_UHetero1sc_sampledSNP_%d.txt"
se2_template <- "PhenotypeSimulated_singleSNP_residSE2_h2%s_UHetero1sc_sampledSNP_%d.txt"

# Output directory
output_dir <- "/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/preWhiten/power_analysis"

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
load(file = "/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/NORMgenomat.samp.Rda")  # Contains 'g'
g <- NORMAL.genomat.samp_mat
cat(sprintf("  Loaded: %d strains × %d SNPs\n", nrow(g), ncol(g)))

# Load kinship matrix
cat("Loading kinship matrix...\n")
fread(file = "/work/jb621/simHaplo/LD_mapping/ars_SS_u7o3_sims_v2/ars_SS_u7o3_sims/StepWiseSuSiE/invpdf_sampling_causal/gemma/grm_inbred/gemmaGRM.full.sXX.txt") -> kinship
kinship_matrix <- as.matrix(kinship)
cat(sprintf("  Loaded: %d × %d\n", nrow(kinship_matrix), ncol(kinship_matrix)))

# Load causal SNP mapping (if available)
cat("Looking for causal SNP metadata...\n")
causal_snp_map <- NULL

metadata_files <- c(
  file.path(base_dir, "causal_snps.txt"),
  file.path(base_dir, "simulation_metadata.txt"),
  "/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/sampled_snp_info.txt"
)

for (meta_file in metadata_files) {
  if (file.exists(meta_file)) {
    cat(sprintf("  Found: %s\n", meta_file))
    causal_snp_map <- read.table(meta_file, header = TRUE, stringsAsFactors = FALSE)
    break
  }
}

if (is.null(causal_snp_map)) {
  cat("  No metadata found - will use SNP column indices\n")
}

# Get strain names
strain_names <- rownames(g)

# ============================================================================
# Verify File Paths and Patterns
# ============================================================================

cat("\n=== Verifying File Patterns ===\n\n")

# Test with first h2 level and first simulation
test_h2 <- h2_levels[1]
test_sim <- 1

# Source the h2_to_int function from power_analysis_functions.R
# source("power_analysis_functions.R")

#test_h2_int <- h2_to_int(test_h2)
test_h2_int <- "0005"

# Construct test file paths
test_pheno <- file.path(base_dir, sprintf(pheno_template, test_h2_int, test_sim))
test_se2 <- file.path(base_dir, sprintf(se2_template, test_h2_int, test_sim))

cat("Expected file paths:\n")
cat(sprintf("  Pheno: %s\n", test_pheno))
cat(sprintf("  SE2:   %s\n", test_se2))

cat("\nFile existence check:\n")
cat(sprintf("  Pheno exists: %s\n", file.exists(test_pheno)))
cat(sprintf("  SE2 exists:   %s\n", file.exists(test_se2)))

if (file.exists(test_pheno)) {
  cat("\nReading test phenotype file...\n")
  test_data <- read.table(test_pheno, header = FALSE, nrows = 5)
  cat(sprintf("  Dimensions: %d rows × %d columns\n", nrow(test_data), ncol(test_data)))
  cat("  First few rows:\n")
  print(test_data)
} else {
  cat("\nWARNING: Test file not found!\n")
  cat("Checking directory contents...\n")
  
  if (dir.exists(base_dir)) {
    all_files <- list.files(base_dir, pattern = "^Phenotype.*\\.txt$")
    cat(sprintf("Found %d phenotype files in directory\n", length(all_files)))
    if (length(all_files) > 0) {
      cat("Sample files (first 5):\n")
      print(head(all_files, 5))
      
      # Try to extract h2 pattern from actual files
      cat("\nTrying to detect h2 pattern from filenames...\n")
      h2_patterns <- str_extract(all_files, "h2[0-9]+")
      unique_h2s <- unique(h2_patterns[!is.na(h2_patterns)])
      cat("Detected h2 codes:", paste(unique_h2s, collapse = ", "), "\n")
    }
  } else {
    cat(sprintf("ERROR: Directory does not exist: %s\n", base_dir))
  }
}

cat("\n")
readline(prompt = "Press [enter] to continue with compilation or Ctrl+C to abort and fix paths...")

# ============================================================================
# Compile Power Simulations (if not already compiled)
# ============================================================================

# Check if already compiled
compiled_file <- file.path(output_dir, "FITMEAN_power_simulations_compiled.Rda")

if (file.exists(compiled_file)) {
  cat("\n=== Loading Pre-compiled Power Simulations ===\n\n")
  load(compiled_file)
  cat(sprintf("Loaded %d observations\n", nrow(power_data)))
  
} else {
  cat("\n=== Compiling Power Simulation Files ===\n\n")
  
  power_data <- compile_power_simulations(
    base_dir = base_dir,
    h2_levels = h2_levels,
    n_sims_per_h2 = n_sims_per_h2,
    file_pattern_pheno = pheno_template,
    file_pattern_se2 = se2_template,
    strain_names = strain_names,
    n_cores = n_cores
  )
  
  # Save compiled data
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  save(power_data, file = compiled_file)
  cat(sprintf("\nSaved compiled data to: %s\n", compiled_file))
}

load(compiled_file)

plotting_data <- power_data %>% 
  filter(trait == "Trait_h2_0025_sim_3")
plot(x = plotting_data$value, y = plotting_data$se2, 
     xlab = "Average Simulated Trait Value per Strain", ylab = "Variance in Simulated Trait Value per Strain",
     main = "Mean-Variance Relationship for Simulated Trait Data with uniformly random variance: h2 = 0.025, sim 1")

# 1. Scatter plot: Mean-Variance Relationship
# p1 <- ggplot(plotting_data, aes(x = value, y = se2)) +
#   geom_point(size = 2) +
#   labs(
#     x     = "Average Simulated Trait Value per Strain",
#     #y     = "Fitted Variance in Simulated Trait Value per Strain",
#     y     = "Variance in Simulated Trait Value per Strain",
#     title = "Mean-Variance Relationship for Simulated Trait Data with uniformly random variance: h2 = 0.025, sim 1"
#   )   +
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

# Scatterplot for the fitted mean-variance relationship
# p4 <- ggplot(plotting_data, aes(x = value, y = se2)) +
#   geom_point() +
#   geom_smooth(
#     method  = "loess",
#     span    = 0.95,
#     color   = "firebrick",
#     fill    = "firebrick",
#     alpha   = 0.15,         # shaded CI ribbon
#     linewidth = 1
#   ) +
#   labs(
#     x     = "Average Simulated Trait Value per Strain",
#     y     = "Fitted Variance in Simulated Trait Value per Strain",
#     title = "Mean-Variance Relationship for Simulated Trait Data with fitted mean-var relationship: h2 = 0.025, sim 1"
#   ) +
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
# p4

# ============================================================================
# Run Power Tests for Both Methods
# ============================================================================

cat("\n=== Running Power Tests ===\n\n")

# Check if results already exist
results_file <- file.path(output_dir, "GEMMA_DefvModTesting_UHetero1sc_NOSHRINKAGE_power_results.Rda")

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
    causal_snp_map = causal_snp_map,
    alpha = alpha,
    n_cores = n_cores
  )
  
  # Save results
  save(power_results, file = results_file)
  cat(sprintf("\nSaved results to: %s\n", results_file))
}

# ============================================================================
# Summarize Power by H2 Level
# ============================================================================

cat("\n=== Summarizing Power by H2 Level ===\n\n")

power_summary <- summarize_power_by_h2(power_results, alpha = alpha)

cat("Power Summary:\n")
print(power_summary, width = Inf)

# ============================================================================
# Calculate Detailed Power Gain Statistics
# ============================================================================

cat("\n=== Calculating Power Gain Statistics ===\n\n")

power_gain_stats <- calculate_power_gain_statistics(power_results)

# Save statistics
stats_file <- file.path(output_dir, "GEMMA_DefVModTesting_NOSHRINKAGE_power_gain_statistics.Rda")
save(power_gain_stats, file = stats_file)

# ============================================================================
# Create Comprehensive Report
# ============================================================================

# Create the mapping
h2_map <- c(
  "0005" = 0.005,
  "001"  = 0.01,
  "0015" = 0.015,
  "002"  = 0.02,
  "0025" = 0.025,
  "003"  = 0.03,
  "0035" = 0.035,
  "004"  = 0.04,
  "0045" = 0.045,
  "005"  = 0.05
  # add all your levels
)

# Apply to power_results
power_results <- power_results %>%
  mutate(h2 = h2_map[as.character(h2)])

# Apply to power_summary  
power_summary <- power_summary %>%
  mutate(h2 = h2_map[as.character(h2)])

# Apply to power_data  
power_data <- power_data %>%
  mutate(h2 = h2_map[as.character(h2)]) %>% 
  select(!h2_int)

# Verify
print(sort(unique(power_results$h2)))

plots <- create_power_report(
  power_results = power_results,
  power_summary = power_summary,
  output_dir = output_dir
)
plots[4]

library(dplyr)
library(tidyr)

power_results_long <- power_results %>%
  # Pivot the two detected columns into long format, creating a 'method' column
  pivot_longer(
    cols = c(detected_default, detected_modified),
    names_to = "method",
    names_prefix = "detected_",
    values_to = "detected"
  ) %>%
  mutate(method = case_when(
    method == "default"  ~ "GEMMA_default",
    method == "modified" ~ "GEMMA_modified"
  )) %>%
  # Each row is one simulation: the causal SNP is a true positive if detected,
  # false negative if not. There are no non-causal SNPs tested here, so FP = 0.
  mutate(
    true_positives  = as.integer(detected),
    false_negatives = as.integer(!detected),
    false_positives = 0L          # no non-causal SNPs in this per-sim structure
  )

summary_ci <- power_results_long %>%
  group_by(method, h2) %>%
  summarise(
    total_TP = sum(true_positives),
    total_FN = sum(false_negatives),
    total_FP = sum(false_positives),
    n_power  = total_TP + total_FN,
    n_fdr    = total_TP + total_FP,
    power_pooled = total_TP / n_power,
    fdr_pooled   = if_else(n_fdr > 0, total_FP / n_fdr, NA_real_),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    power_ci    = list(binom::binom.confint(total_TP, n_power, method = "wilson")),
    power_lower = power_ci$lower,
    power_upper = power_ci$upper,
    fdr_ci      = list(binom::binom.confint(total_FP, max(n_fdr, 1L), method = "wilson")),
    fdr_lower   = fdr_ci$lower,
    fdr_upper   = fdr_ci$upper
  )

p1 <- ggplot(summary_ci, aes(x = h2, y = power_pooled, color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = power_lower, ymax = power_upper), width = 0.002) +
  geom_hline(yintercept = 0.8, linetype = "dashed", 
             color = "gray50", alpha = 0.5) +
  annotate("text", x = min(summary_ci$h2), y = 0.8, 
           label = "80% Power", vjust = -0.5, hjust = -0.25, color = "gray50", size = 7) +
  scale_color_discrete(labels = c("GEMMA_default"  = "GEMMA",
                                  "GEMMA_modified" = "GEMMAmod")) +
  labs(
    title    = "Power Comparison: GEMMA vs GEMMAmod",
    #subtitle = "Uniform Variance",
    subtitle = "Fitted Mean-Variance",
    x        = "Heritability (h²)",          # <-- updated from Fold Change
    y        = "Power (Correct Causal SNP Detection)",
    color    = "Method"
  ) +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(size = 14),
    axis.text.y      = element_text(size = 14),
    plot.title       = element_text(hjust = 0.5, size = 22),
    plot.subtitle    = element_text(hjust = 0.5, size = 14),
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 14),
    axis.title.x     = element_text(size = 18),
    axis.title.y     = element_text(size = 18)
  )
p1

plots_file <- file.path(output_dir, "GEMMA_DefvModTesting_UHetero1sc_NOSHRINKAGE_power_plots.Rda")

save(plots, file = plots_file)

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

