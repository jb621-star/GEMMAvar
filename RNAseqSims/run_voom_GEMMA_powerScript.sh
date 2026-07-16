#!/bin/bash
#SBATCH --job-name=gemma_eqtl_sim
#SBATCH --output=gemma_eqtl_sim_%j.out
#SBATCH --error=gemma_eqtl_sim_%j.err
#SBATCH --time=12-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --partition=scavenger

# Load required modules (adjust these based on your cluster)
module load R/4.4.3

# Set number of threads for parallel processing
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Start time: $(date)"
echo "Working directory: $SLURM_SUBMIT_DIR"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory allocated: $SLURM_MEM_PER_NODE MB"

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Run the R script
Rscript /work/jb621/GEMMA_2_pickrell_sims/limma_voom_to_GEMMA_eQTL_powerScript_wCov_fc2_03302026.R

# Print completion information
echo "End time: $(date)"
echo "Job completed successfully"
