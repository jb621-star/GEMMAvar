#!/bin/bash
## ============================================================
## 03_submit_gemma_array.sh
## SLURM array job — one task per peak
##
## Submit with:
##   sbatch 03_submit_gemma_array.sh
##
## To limit concurrent jobs (e.g. 500 at a time):
##   sbatch --array=1-$(wc -l < gemma_peak_inputs/peak_list.txt)%500 03_submit_gemma_array.sh
## ============================================================

#SBATCH --account=pro00117530
#SBATCH --job-name=gemmamod_caqtl
#SBATCH --output=logs/gemmamod_caqtl_%A_%a.out
#SBATCH --error=logs/gemmamod_caqtl_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

# array-size needs to be specified within the sbatch command

module load miniforge
conda activate gemma_r422

# Create log directory if it doesn't exist
mkdir -p logs

# Apply offset if set (for batches beyond MaxArraySize)
OFFSET=${OFFSET:-0}
ACTUAL_IDX=$((SLURM_ARRAY_TASK_ID + OFFSET))
 
# Get peak ID for this array task
PEAK=$(sed -n "${ACTUAL_IDX}p" caqtl_gene_inputs/peak_list_rerun.txt)
 
if [ -z "$PEAK" ]; then
  echo "No peak found for index ${ACTUAL_IDX} (task ${SLURM_ARRAY_TASK_ID} + offset ${OFFSET}), exiting."
  exit 1
fi
 
echo "Task ${SLURM_ARRAY_TASK_ID} | Offset ${OFFSET} | Index ${ACTUAL_IDX} | peak: ${PEAK}"
 
Rscript run_peak_gemmamod_02.R "${PEAK}"
 
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "ERROR: Rscript failed for peak ${PEAK} with exit code ${EXIT_CODE}"
  exit $EXIT_CODE
fi
 
echo "Completed: ${PEAK}"