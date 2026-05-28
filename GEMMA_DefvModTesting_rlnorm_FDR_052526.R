library(Matrix)
library(tidyverse)
library(data.table)
library(parallel)
library(data.table)

#load("/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/normal_null_m25_sd5_1000sims.Rda")
#load("/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/normal_null_m25_within3_between4_1000sims.Rda")

# Base directory
base_dir <- "/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/preWhiten"

# Number of simulations
n_sims <- 1000

# File path templates
pheno_template <- file.path(base_dir, "PhenotypeSimulated_singleSNP_null_rlnormHETERO_mneg0sd15_0527_sampledSNP_%d.txt")
se2_template <- file.path(base_dir, "PhenotypeSimulated_singleSNP_residSE2_null_rlnormHETERO_mneg0sd15_0527_sampledSNP_%d.txt")

# pheno_template <- file.path(base_dir, "PhenotypeSimulated_singleSNP_null_FITmeanvar_sampledSNP_%d.txt")
# se2_template <- file.path(base_dir, "PhenotypeSimulated_singleSNP_residSE2_null_FITmeanvar_sampledSNP_%d.txt")

#pheno_template <- file.path(base_dir, "PhenotypeSimulated_singleSNP_null_UHetero1scOUTLIERS_0527_sampledSNP_%d.txt")
#se2_template <- file.path(base_dir, "PhenotypeSimulated_singleSNP_residSE2_null_UHetero1scOUTLIERS_0527_sampledSNP_%d.txt")


# Inputs
load(file = "/work/jb621/simHaplo/LD_mapping/taggingSNPs_v4_geno_df.Rda")
#kinship <- read.table("/work/jb621/simHaplo/LD_mapping/ars_SS_u7o3_sims_v2/ars_SS_u7o3_sims/StepWiseSuSiE/invpdf_sampling_causal/gemma/grm_inbred/gemmaGRM.full.sXX.txt", sep="\t", header=FALSE)
kinship_matrix <- as.matrix(fread("/work/jb621/simHaplo/LD_mapping/ars_SS_u7o3_sims_v2/ars_SS_u7o3_sims/StepWiseSuSiE/invpdf_sampling_causal/gemma/grm_inbred/gemmaGRM.full.sXX.txt"))
#cat("Dimensions:    ", dim(kinship_matrix), "\n")   # should be 192 192
#cat("Diagonal range:", range(diag(kinship_matrix)), "\n")  # should be numeric ~0.3-5
#cat("Is numeric:    ", is.numeric(kinship_matrix), "\n")   # must be TRUE

load(file = "/hpc/group/baughlab/jb621/GEMMA_2_testing/sampledSNPs/NORMgenomat.samp.Rda")

#fread(file = "/hpc/group/baughlab/jb621/GEMMA_2_testing/slope100/output/gemmaGRM.NORMAL.sXX.txt") -> kinship2
#load(file = "/work/jb621/simHaplo/LD_mapping/taggingSNPs_v4_geno_df.Rda")
# ============================================================================
# Power Analysis for GEMMA Methods
# ============================================================================
# Functions to calculate and compare power between Default and Modified GEMMA
# across different heritability levels
# ============================================================================
# 1. Check what kinship_matrix actually looks like
cat("Dimensions:", dim(kinship_matrix), "\n")
cat("Class:", class(kinship_matrix), "\n")
cat("Is numeric:", is.numeric(kinship_matrix), "\n")
cat("Diagonal range:", range(diag(kinship_matrix)), "\n")
cat("Off-diagonal range:", range(kinship_matrix[lower.tri(kinship_matrix)]), "\n")
print(kinship_matrix[1:5, 1:5])

# 2. Check what the phenotype and genotype look like for one sim
s <- 1
sim_data  <- simulated_null_data[simulated_null_data$sim_id == s, ]
causal_snp <- metadata$causal_snp[metadata$sim_id == s]
phenotype  <- sim_data$value
snp_geno   <- as.matrix(NORMAL.genomat.samp_mat[, causal_snp, drop = FALSE])

cat("\nPhenotype length:", length(phenotype), "\n")
cat("Phenotype variance:", var(phenotype), "\n")
cat("SNP genotype unique values:", unique(snp_geno), "\n")
cat("SNP genotype variance:", var(snp_geno), "\n")

# 3. Run one simulation manually and inspect intermediate objects
p1 <- run_default_gemma(phenotype, snp_geno, kinship_matrix)
cat("\nSingle p-value (default):", p1, "\n")

# 4. Check p-value distribution -- is it uniform or something else?
cat("\nQuantiles of p_default:\n")
print(quantile(results_df$p_default, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))

# 5. Are there NAs?
cat("\nNA count - default:", sum(is.na(results_df$p_default)), "\n")
cat("NA count - modified:", sum(is.na(results_df$p_modified)), "\n")

# 6. Check if phenotypes were simulated with kinship structure
cat("\nCorrelation between phenotype and first kinship PC:\n")
eig <- eigen(kinship_matrix, symmetric = TRUE)
pc1 <- eig$vectors[, 1]
print(cor(phenotype, pc1))

# Distribution of per-SNP correlations with phenotype across all 1000 sims
snp_pheno_cors <- sapply(1:1000, function(s) {
  sim_data   <- simulated_null_data[simulated_null_data$sim_id == s, ]
  causal_snp <- metadata$causal_snp[metadata$sim_id == s]
  phenotype  <- sim_data$value
  snp_geno   <- as.numeric(NORMAL.genomat.samp_mat[, causal_snp])
  cor(phenotype, snp_geno)
})

hist(snp_pheno_cors, breaks = 40,
     main = "SNP-phenotype correlations across 1000 null sims",
     xlab = "Pearson r")

# Check correlation between p-value and SNP-phenotype correlation
plot(abs(snp_pheno_cors), -log10(results_df$p_default),
     pch = 20, xlab = "|r(SNP, phenotype)|",
     ylab = "-log10(p)", main = "Is inflation driven by SNP-phenotype correlation?")

# Re-run a sample of sims but extract lambda too
lambda_check <- sapply(1:100, function(s) {
  sim_data   <- simulated_null_data[simulated_null_data$sim_id == s, ]
  causal_snp <- metadata$causal_snp[metadata$sim_id == s]
  phenotype  <- as.matrix(sim_data$value)
  snp_geno   <- as.matrix(NORMAL.genomat.samp_mat[, causal_snp, drop = FALSE])
  W          <- matrix(1, nrow = length(phenotype), ncol = 1)
  
  data <- prepare_gemma_data(pheno = phenotype, geno = snp_geno,
                             covars = W, kinship = kinship_matrix)
  
  eig     <- eigen(data$K, symmetric = TRUE)
  evalues <- pmax(eig$values, 0)
  Yr <- t(eig$vectors) %*% data$Y
  Wr <- t(eig$vectors) %*% data$W
  Xr <- t(eig$vectors) %*% data$X
  
  lambda <- calc_lambda_restricted(evalues, Yr, cbind(Wr, Xr), grid = TRUE)
  return(lambda)
})

cat("Lambda range:", range(lambda_check), "\n")
cat("Lambda < 1e-3:", sum(lambda_check < 1e-3), "\n")
cat("Lambda > 1e3: ", sum(lambda_check > 1e3), "\n")
hist(log10(lambda_check), breaks = 30,
     main = "Distribution of log10(lambda)",
     xlab = "log10(lambda)")

# Cross-reference with p-values
plot(log10(lambda_check), -log10(results_df$p_default[1:100]),
     pch = 20, xlab = "log10(lambda)",
     ylab = "-log10(p)",
     main = "Does lambda predict significance?")

lambda_recheck <- sapply(1:10, function(s) {
  sim_data   <- simulated_null_data[simulated_null_data$sim_id == s, ]
  causal_snp <- metadata$causal_snp[metadata$sim_id == s]
  phenotype  <- as.matrix(sim_data$value)
  snp_geno   <- as.matrix(NORMAL.genomat.samp_mat[, causal_snp, drop = FALSE])
  W          <- matrix(1, nrow = length(phenotype), ncol = 1)
  
  data    <- prepare_gemma_data(pheno = phenotype, geno = snp_geno,
                                covars = W, kinship = kinship_matrix)
  eig     <- eigen(data$K, symmetric = TRUE)
  evalues <- pmax(eig$values, 0)
  Yr <- t(eig$vectors) %*% data$Y
  Wr <- t(eig$vectors) %*% data$W
  Xr <- t(eig$vectors) %*% data$X
  
  calc_lambda_restricted(evalues, Yr, cbind(Wr, Xr), grid = TRUE)
})

cat("Lambda range after fix:", range(lambda_recheck), "\n")
cat("Lambda < 1e-3:         ", sum(lambda_recheck < 1e-3), "\n")

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
  select(!causal_snp) %>%       # drop the NA placeholder added by compile_simulated_data
  left_join(metadata, by = "sim_id")  # attach the real causal SNP labels

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
      grid = TRUE,
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
      grid = TRUE,
      nproc = 1
    )
  }, error = function(e) {
    return(data.frame(p_wald = NA_real_))
  })
  
  return(result$p_wald[1])
}

results_list <- list()
#test_kinship <- diag(1, nrow = n, ncol = n)
for (s in 1:1000) {
  #s <- 1
  # ---- Subset this simulation's data ----
  sim_data <- simulated_null_data[simulated_null_data$sim_id == s, ]
  
  causal_snp  <- metadata$causal_snp[metadata$sim_id == s]
  phenotype   <- sim_data$value
  pheno_var   <- sim_data$se2
  snp_geno    <- as.matrix(NORMAL.genomat.samp_mat[, causal_snp, drop = FALSE])
  
  # ---- Run both GEMMA variants ----
  p_default  <- tryCatch(
    run_default_gemma(phenotype, snp_geno, kinship_matrix),
    error = function(e) NA_real_
  )
  
  p_modified <- tryCatch(
    run_modified_gemma(phenotype, snp_geno, kinship_matrix,
                       phenotype_variance = pheno_var,
                       lambda_shrink = 0),
    error = function(e) NA_real_
  )
  
  results <- data.frame(
    sim_id     = s,
    p_default  = p_default,
    p_modified = p_modified
  )
  
  results_list[[s]] <- results 
}

results_df <- bind_rows(results_list, .id = "sim_id")
results_df %>% 
  summarise(def_FPR = sum(p_default < 0.05) / 1000,
            mod_FPR = sum(p_modified < 0.05) / 1000)

table(results_df$p_default < 0.05)


# ============== Trying to diagnose what was wrong with the simulated phenotypes
# ============== CONCLUSION: IT'S THE KINSHIOP (?)...why wasn't this the issue yesterday????
random_snps <- sample(colnames(NORMAL.genomat.samp_mat))

perm_results_list <- list()
for (s in 1:1000) {
  
  sim_data <- simulated_null_data[
    simulated_null_data$sim_id == s, ]
  
  phenotype <- sim_data$value
  
  snp_geno <- as.matrix(
    NORMAL.genomat.samp_mat[
      , random_snps[s], drop = FALSE]
  )
  
  p_default <- run_default_gemma(
    phenotype,
    snp_geno,
    kinship_matrix
  )
  
  perm_results_list[[s]] <- p_default
}

results_df <-unlist(perm_results_list)
sum(results_df < 0.05) / 1000

run_default_gemma <- function(
    phenotype,
    snp_genotype,
    kinship
) {
  
  data <- prepare_gemma_data(
    pheno = phenotype,
    geno = snp_genotype,
    covars = NULL,
    kinship = kinship
  )
  
  result <- pygemma(
    Y = data$Y,
    X = data$X,
    W = data$W,
    K = data$K,
    verbose = 0,
    grid = FALSE,
    nproc = 1
  )
  
  result$p_wald[1]
}

# Purely synthetic genotype test
set.seed(1)

y_test <- rnorm(192)

pvals <- replicate(1000, {
  
  rand_snp <- sample(
    colnames(NORMAL.genomat.samp_mat), 1
  )
  
  snp_geno <- as.matrix(
    NORMAL.genomat.samp_mat[, rand_snp, drop = FALSE]
  )
  
  snp_geno <- matrix(
    as.numeric(snp_geno),
    ncol = 1
  )
  
  run_default_gemma(
    y_test,
    snp_geno,
    kinship_matrix
  )
})

sum(pvals < 0.05, na.rm = TRUE)

run_lm_only <- function(y, x) {
  
  data <- prepare_gemma_data(
    pheno = y,
    geno = x,
    covars = NULL,
    kinship = diag(length(y))
  )
  
  result <- pygemma(
    Y = data$Y,
    X = data$X,
    W = data$W,
    K = diag(length(y)),
    verbose = 0,
    grid = FALSE,
    nproc = 1
  )
  
  result$p_wald[1]
}

set.seed(1)

pvals <- replicate(1000, {
  
  rand_snp <- sample(
    colnames(NORMAL.genomat.samp_mat), 1
  )
  
  x <- matrix(
    as.numeric(
      NORMAL.genomat.samp_mat[, rand_snp]
    ),
    ncol = 1
  )
  
  run_lm_only(y_test, x)
})

sum(pvals < 0.05)

# lambda shrinkage turned the null FDR from 0.024 to 0.055
# find out what exactly this means for everything else
hist(results_df$p_default)
hist(results_df$p_modified)

# Build 2x2 contingency table of rejection patterns
results_df <- results_df %>%
  mutate(
    rej_default  = p_default  < 0.05,
    rej_modified = p_modified < 0.05
  )

# Paired contingency table
#                  mod_rejected  mod_not_rejected
# def_rejected         a               b
# def_not_rejected     c               d
# Only b and c matter for McNemar (the discordant pairs)

mc_table <- table(
  default  = results_df$rej_default,
  modified = results_df$rej_modified
)
print(mc_table)

mcnemar.test(mc_table)

# Are each individually uniform?
ks_default  <- ks.test(results_df$p_default,  "punif")
ks_modified <- ks.test(results_df$p_modified, "punif")

cat("KS test vs. uniform:\n")
cat(sprintf("  Default GEMMA:  D = %.4f, p = %.4f\n", ks_default$statistic,  ks_default$p.value))
cat(sprintf("  Modified GEMMA: D = %.4f, p = %.4f\n", ks_modified$statistic, ks_modified$p.value))

library(ggplot2)

results_df_long <- results_df %>%
  pivot_longer(cols = c(p_default, p_modified),
               names_to = "method",
               values_to = "p_value") %>%
  mutate(method = case_when(method == "p_default"  ~ "Default GEMMA",
                            method == "p_modified"  ~ "Modified GEMMA"))

ggplot(results_df_long, aes(x = p_value, fill = method, color = method)) +
  geom_histogram(aes(y = after_stat(count)),
                 bins = 20, alpha = 0.5, position = "identity",
                 linewidth = 0.3) +
  geom_hline(yintercept = 1000 / 20,
             linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_fill_manual(values  = c("Default GEMMA" = "#378ADD",
                                "Modified GEMMA" = "#D85A30")) +
  scale_color_manual(values = c("Default GEMMA" = "#378ADD",
                                "Modified GEMMA" = "#D85A30")) +
  annotate("text", x = 0.85, y = (1000 / 20) + 3,
           label = "Uniform null", size = 3.5, color = "black") +
  labs(x = "p-value", y = "Count", fill = NULL, color = NULL) +
  theme_classic() +
  theme(legend.position = "top")

# 1. Check the realized vs. stored variance
#    Do the se2 values actually match per-individual residual variance?
realized_var <- sapply(1:1000, function(s) {
  sim_data <- simulated_null_data[simulated_null_data$sim_id == s, ]
  var(sim_data$value)  # should be ~1 after scaling
})
hist(realized_var)  # should be tightly around 1

# 2. Check the weight distribution going into GEMMAmod
#    Are weights sensible or extreme?
check_weights <- function(se2_vec, lambda = 0.2) {
  log_s2 <- log(se2_vec)
  log_s2_shr <- (1 - lambda) * log_s2 + lambda * mean(log_s2)
  w <- 1 / sqrt(exp(log_s2_shr))
  cat("Weight range:", round(range(w), 3), "\n")
  cat("Weight CV:", round(sd(w)/mean(w), 3), "\n")  # >0.5 suggests strong weighting
}

sim1 <- simulated_null_data[simulated_null_data$sim_id == 1, ]
check_weights(sim1$se2)

# 3. PP plot to distinguish conservatism from miscalibration
plot(sort(results_df$p_default), 
     ppoints(1000), pch = 20, col = "blue",
     xlab = "Observed p", ylab = "Expected p (uniform)",
     main = "PP plot: null calibration")
points(sort(results_df$p_modified), 
       ppoints(1000), pch = 20, col = "red")
abline(0, 1, lty = 2)
legend("topleft", c("Default", "Modified"), col = c("blue","red"), pch = 20)

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
             label = "Equal Power", vjust = -0.5, color = "red", size = 4) +
    scale_x_continuous(breaks = unique(power_summary$h2)) +
    labs(title = "Power Ratio: Modified / Default GEMMA",
         subtitle = "Values > 1 indicate Modified has higher power",
         x = expression(Heritability~(h^2)),
         y = "Power Ratio (Modified / Default)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  
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
                  width = 0.02, size = 0.8) +
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
# Option 1: Compile Power Simulations (if not already compiled)
# ============================================================================

# Check if already compiled
compiled_file <- file.path(output_dir, "power_simulations_compiled.Rda")

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

plotting_data <- power_data %>% 
  filter(trait == "Trait_h2_0025_sim_1")
plot(x = plotting_data$value, y = plotting_data$se2, 
     xlab = "Average Simulated Trait Value per Strain", ylab = "Variance in Simulated Trait Value per Strain",
     main = "Mean-Variance Relationship for Simulated Trait Data with uniformly random variance: h2 = 0.025, sim 1")

# ============================================================================
# Run Power Tests for Both Methods
# ============================================================================

cat("\n=== Running Power Tests ===\n\n")

# Check if results already exist
results_file <- file.path(output_dir, "GEMMA_DefvModTesting_UHetero1sc_power_results.Rda")

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
stats_file <- file.path(output_dir, "GEMMA_DefVModTesting_power_gain_statistics.Rda")
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

power_results_uhetero <- power_results

# Apply to power_summary  
power_summary <- power_summary %>%
  mutate(h2 = h2_map[as.character(h2)])

power_summary_uhetero <- power_summary

# Apply to power_data  
power_data <- power_data %>%
  mutate(h2 = h2_map[as.character(h2)]) %>% 
  select(!h2_int)

power_data_uhetero <- power_data

# Verify
print(sort(unique(power_results$h2)))

plots <- create_power_report(
  power_results = power_results,
  power_summary = power_summary,
  output_dir = output_dir
)
plots[1]

plots_file <- file.path(output_dir, "GEMMA_DefvModTesting_UHetero1sc_power_plots.Rda")

save(plots, file = plots_file)
# ============================================================================
# Additional Analyses
# ============================================================================

cat("\n=== Additional Analyses ===\n\n")

# 1. P-value comparison scatter plots by h2
cat("Creating p-value comparison plots...\n")

for (h2_val in unique(power_results$h2)) {
  h2_data <- power_results %>% filter(h2 == h2_val)
  
  p <- ggplot(h2_data, aes(x = log10p_default, y = log10p_modified)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = -log10(alpha), linetype = "dashed", alpha = 0.5) +
    labs(title = sprintf("P-value Comparison at h² = %.2f", h2_val),
         x = "-log10(p) Default GEMMA",
         y = "-log10(p) Modified GEMMA") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  p
  
  # ggsave(file.path(output_dir, sprintf("pvalue_scatter_h2_%.2f.png", h2_val)),
  #        p, width = 8, height = 8, dpi = 300)
}
cat(sprintf("  Saved %d scatter plots\n", length(unique(power_results$h2))))

library(tidyverse)

# ── Step 1: Build a SNP name lookup from the genotype matrix ──────────────────
# Assumes colnames(g) are in the same "I:62576" format as power_results$causal_snp
snp_names <- colnames(g)       # character vector of length 1000
strain_names <- rownames(g)    # character vector of length 192

# ── Step 2: For each simulation, extract the causal SNP genotype per strain ───
# power_results has one row per simulation — join causal SNP allele onto power_data

causal_alleles <- power_results %>%
  select(sim_id, trait, causal_snp, h2, detected_default, detected_modified) %>%
  rowwise() %>%
  mutate(
    snp_col = which(snp_names == causal_snp),  # index into g
    allele_vec = list(
      tibble(
        strain = strain_names,
        allele = g[, snp_col]  # 0 = REF, 1 = ALT
      )
    )
  ) %>%
  ungroup() %>%
  unnest(allele_vec)

# ── Step 3: Join onto power_data ──────────────────────────────────────────────
# power_data has: strain, value (strain mean), se2 (within-strain var), sim_id, trait, h2

plot_data <- power_data %>%
  left_join(causal_alleles, by = c("sim_id", "trait", "h2", "strain"))

# ── Step 4: Compute per-simulation summary stats ──────────────────────────────
sim_summaries <- plot_data %>%
  group_by(sim_id, trait, h2, detected_default, detected_modified, allele) %>%
  summarise(
    group_mean    = mean(value, na.rm = TRUE),
    group_mean_se2 = mean(se2, na.rm = TRUE),
    n_strains     = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = allele,
    values_from = c(group_mean, group_mean_se2, n_strains),
    names_sep   = "_allele"
  ) %>%
  mutate(
    mean_diff = group_mean_allele2 - group_mean_allele0,       # ALT - REF mean
    var_ratio = group_mean_se2_allele2 / group_mean_se2_allele0  # ALT / REF within-strain var
  )

# Optional: join back skewness/kurtosis computed from strain means per simulation
skew_stats <- plot_data %>%
  group_by(sim_id, trait, h2) %>%
  summarise(
    skewness = moments::skewness(value, na.rm = TRUE),
    kurtosis = moments::kurtosis(value, na.rm = TRUE),
    maf      = { a <- mean(allele, na.rm = TRUE); min(a, 1 - a) },
    .groups  = "drop"
  )

sim_summaries <- sim_summaries %>%
  left_join(skew_stats, by = c("sim_id", "trait", "h2"))

# ── Step 5: Plots ──────────────────────────────────────────────────────────────

# -- P1: Core detectability landscape (use detected_default or detected_modified)
p1 <- ggplot(sim_summaries, aes(x = mean_diff, y = var_ratio, color = detected_modified)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c("TRUE" = "#2166ac", "FALSE" = "#d73027"),
    labels = c("TRUE" = "Detected", "FALSE" = "Not detected")
  ) +
  facet_wrap(~ h2, labeller = label_both) +
  labs(
    x     = "Mean difference between ALT and REF strains",
    y     = "Within-strain variance ratio (ALT / REF)",
    color = "QTL detected",
    title = "Detectability landscape across 1,000 simulations"
  ) +
  theme_bw()
p1

# -- P2: Facet by skewness quartile to see distributional effects
sim_summaries <- sim_summaries %>%
  mutate(skew_quartile = ntile(abs(skewness), 4) %>%
           factor(labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high skew)")))

p2 <- ggplot(sim_summaries, aes(x = mean_diff, y = var_ratio, color = detected_modified)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c("TRUE" = "#2166ac", "FALSE" = "#d73027"),
    labels = c("TRUE" = "Detected", "FALSE" = "Not detected")
  ) +
  facet_wrap(~ skew_quartile, nrow = 1) +
  labs(
    x     = "Mean difference (ALT - REF)",
    y     = "Within-strain variance ratio (ALT / REF)",
    color = "QTL detected",
    title = "Detectability by skewness quartile of strain means"
  ) +
  theme_bw()
p2

# -- P3: Detection rate binned by mean_diff, colored by method
sim_summaries %>%
  pivot_longer(cols = c(detected_default, detected_modified),
               names_to = "method", values_to = "detected") %>%
  mutate(mean_diff_bin = cut(mean_diff, breaks = 15)) %>%
  group_by(mean_diff_bin, method) %>%
  summarise(detection_rate = mean(detected), n = n(), .groups = "drop") %>%
  ggplot(aes(x = mean_diff_bin, y = detection_rate, fill = method)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("detected_default" = "#4dac26", "detected_modified" = "#2166ac")) +
  labs(
    x = "Mean difference ALT - REF (binned)",
    y = "Detection rate",
    title = "Power as a function of allele group mean difference"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -- P4: Detection rate binned by variance ratio
sim_summaries %>%
  pivot_longer(cols = c(detected_default, detected_modified),
               names_to = "method", values_to = "detected") %>%
  mutate(var_ratio_bin = cut(var_ratio, breaks = 15)) %>%
  group_by(var_ratio_bin, method) %>%
  summarise(detection_rate = mean(detected), n = n(), .groups = "drop") %>%
  ggplot(aes(x = var_ratio_bin, y = detection_rate, fill = method)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("detected_default" = "#4dac26", "detected_modified" = "#2166ac")) +
  labs(
    x = "Within-strain variance ratio ALT/REF (binned)",
    y = "Detection rate",
    title = "Power as a function of within-strain variance heterogeneity: Uniform"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================================
# Create Summary Report Document
# ============================================================================

cat("\n=== Creating Summary Report ===\n\n")

report_file <- file.path(output_dir, "POWER_ANALYSIS_REPORT.txt")

sink(report_file)
cat("=" , rep("=", 78), "=\n", sep = "")
cat("POWER ANALYSIS REPORT: Default vs Modified GEMMA\n")
cat("=" , rep("=", 78), "=\n\n", sep = "")

cat(sprintf("Analysis Date: %s\n", Sys.time()))
cat(sprintf("Significance Threshold: α = %.3f\n", alpha))
cat(sprintf("H² Levels Tested: %s\n", paste(h2_levels, collapse = ", ")))
cat(sprintf("Simulations per H²: %d\n", n_sims_per_h2))
cat(sprintf("Total Tests: %d\n\n", nrow(power_results)))

cat("=" , rep("=", 78), "=\n", sep = "")
cat("OVERALL POWER SUMMARY\n")
cat("=" , rep("=", 78), "=\n\n", sep = "")

cat(sprintf("Default GEMMA Power:  %.2f%%\n", 
            100 * mean(power_results$detected_default, na.rm = TRUE)))
cat(sprintf("Modified GEMMA Power: %.2f%%\n", 
            100 * mean(power_results$detected_modified, na.rm = TRUE)))
cat(sprintf("Absolute Gain:        %.2f%%\n",
            100 * (mean(power_results$detected_modified, na.rm = TRUE) - 
                     mean(power_results$detected_default, na.rm = TRUE))))
cat(sprintf("Relative Gain:        %.1f%%\n\n",
            100 * (mean(power_results$detected_modified, na.rm = TRUE) / 
                     mean(power_results$detected_default, na.rm = TRUE) - 1)))

cat("=" , rep("=", 78), "=\n", sep = "")
cat("POWER BY HERITABILITY LEVEL\n")
cat("=" , rep("=", 78), "=\n\n", sep = "")

print(power_summary %>% 
        select(h2, n_tests, power_default, power_modified, power_gain = power_difference), 
      width = Inf)

cat("\n")
cat("=" , rep("=", 78), "=\n", sep = "")
cat("P-VALUE SIGNIFICANCE GAIN\n")
cat("=" , rep("=", 78), "=\n\n", sep = "")

pval_gain <- power_results %>%
  filter(detected_default & detected_modified) %>%
  summarise(
    n = n(),
    mean_p_fold = mean(p_default / p_modified, na.rm = TRUE),
    median_p_fold = median(p_default / p_modified, na.rm = TRUE),
    mean_log10p_gain = mean(log10p_modified - log10p_default, na.rm = TRUE),
    median_log10p_gain = median(log10p_modified - log10p_default, na.rm = TRUE),
    pct_modified_stronger = 100 * mean(log10p_modified > log10p_default, na.rm = TRUE)
  )

cat("Among signals detected by BOTH methods:\n\n")
cat(sprintf("  N = %d\n", pval_gain$n))
cat(sprintf("  Mean p-value fold change:   %.2fx (Modified more significant)\n", 
            pval_gain$mean_p_fold))
cat(sprintf("  Median p-value fold change: %.2fx\n", pval_gain$median_p_fold))
cat(sprintf("  Mean -log10(p) gain:        %.2f\n", pval_gain$mean_log10p_gain))
cat(sprintf("  Median -log10(p) gain:      %.2f\n", pval_gain$median_log10p_gain))
cat(sprintf("  Modified stronger:          %.1f%% of cases\n\n", 
            pval_gain$pct_modified_stronger))

cat("By Heritability Level:\n\n")

h2_pval_gain <- power_results %>%
  filter(detected_default & detected_modified) %>%
  group_by(h2) %>%
  summarise(
    n = n(),
    mean_log10p_gain = mean(log10p_modified - log10p_default, na.rm = TRUE),
    median_log10p_gain = median(log10p_modified - log10p_default, na.rm = TRUE),
    .groups = 'drop'
  )

print(h2_pval_gain, width = Inf)

cat("\n")
cat("=" , rep("=", 78), "=\n", sep = "")
cat("DETECTION PATTERNS\n")
cat("=" , rep("=", 78), "=\n\n", sep = "")

detection_summary <- power_results %>%
  summarise(
    both = sum(detected_default & detected_modified, na.rm = TRUE),
    default_only = sum(detected_default & !detected_modified, na.rm = TRUE),
    modified_only = sum(!detected_default & detected_modified, na.rm = TRUE),
    neither = sum(!detected_default & !detected_modified, na.rm = TRUE),
    total = n()
  ) %>%
  mutate(
    pct_both = 100 * both / total,
    pct_default_only = 100 * default_only / total,
    pct_modified_only = 100 * modified_only / total,
    pct_neither = 100 * neither / total
  )

cat(sprintf("Both methods detected:     %5d (%.1f%%)\n", 
            detection_summary$both, detection_summary$pct_both))
cat(sprintf("Default only detected:     %5d (%.1f%%)\n", 
            detection_summary$default_only, detection_summary$pct_default_only))
cat(sprintf("Modified only detected:    %5d (%.1f%%)\n", 
            detection_summary$modified_only, detection_summary$pct_modified_only))
cat(sprintf("Neither detected:          %5d (%.1f%%)\n\n", 
            detection_summary$neither, detection_summary$pct_neither))

cat("=" , rep("=", 78), "=\n", sep = "")
cat("OUTPUT FILES\n")
cat("=" , rep("=", 78), "=\n\n", sep = "")

cat("Data Files:\n")
cat(sprintf("  - %s\n", "power_results.Rda"))
cat(sprintf("  - %s\n", "power_summary_by_h2.csv"))
cat(sprintf("  - %s\n", "power_results_full.csv"))
cat(sprintf("  - %s\n", "power_gain_statistics.Rda"))

cat("\nPlots:\n")
cat(sprintf("  - %s\n", "power_curves.png"))
cat(sprintf("  - %s\n", "power_ratio.png"))
cat(sprintf("  - %s\n", "power_gain.png"))
cat(sprintf("  - %s\n", "pvalue_significance_gain.png"))
cat(sprintf("  - %s\n", "detection_patterns.png"))
cat(sprintf("  - %s\n", "pvalue_ratio_distribution.png"))
cat(sprintf("  - %s\n", "detection_concordance_heatmap.png"))
cat(sprintf("  - pvalue_scatter_h2_*.png (one per h² level)\n"))

cat("\n")
cat("=" , rep("=", 78), "=\n", sep = "")
cat("END OF REPORT\n")
cat("=" , rep("=", 78), "=\n", sep = "")

sink()

cat(sprintf("Summary report saved to: %s\n", report_file))

# ============================================================================
# Final Summary
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep = "")
cat("POWER ANALYSIS COMPLETE!\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Key Findings:\n")
cat(sprintf("  Overall Power Gain: %.2f%%\n",
            100 * (mean(power_results$detected_modified, na.rm = TRUE) - 
                     mean(power_results$detected_default, na.rm = TRUE))))
cat(sprintf("  Modified %.1f%% more powerful than Default\n",
            100 * (mean(power_results$detected_modified, na.rm = TRUE) / 
                     mean(power_results$detected_default, na.rm = TRUE) - 1)))
cat(sprintf("  P-values %.2fx more significant (median)\n",
            pval_gain$median_p_fold))

cat(sprintf("\nAll outputs saved to: %s\n", output_dir))
cat(sprintf("\nRead full report: %s\n", report_file))

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
calc_lambda_restricted <- function(eigenVals, Y, W, grid = TRUE) {
  
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

