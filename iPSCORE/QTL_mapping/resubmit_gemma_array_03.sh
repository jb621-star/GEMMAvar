#!/bin/bash
## ============================================================
## 03_submit_gemma_array.sh
## SLURM array job — one task per gene
##
## Submit with:
##   sbatch 03_submit_gemma_array.sh
##
## To limit concurrent jobs (e.g. 500 at a time):
##   sbatch --array=1-$(wc -l < gemma_gene_inputs/gene_list.txt)%500 submit_gemma_array_03.sh
## ============================================================

#SBATCH --account=pro00117530
#SBATCH --job-name=gemmamod_eqtl
#SBATCH --output=logs/gemmamod_eqtl_%A_%a.out
#SBATCH --error=logs/gemmamod_eqtl_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00

# array-size needs to be specified within the sbatch command

module load miniforge
conda activate gemma_r422

# Create log directory if it doesn't exist
mkdir -p logs

# Use custom gene list if provided, otherwise default
LIST=${GENE_LIST:-gemma_gene_inputs/gene_list.txt}

ACTUAL_IDX=$((SLURM_ARRAY_TASK_ID + ${OFFSET:-0}))
GENE=$(sed -n "${ACTUAL_IDX}p" "$LIST")
 
if [ -z "$GENE" ]; then
  echo "No gene found for index ${ACTUAL_IDX} (task ${SLURM_ARRAY_TASK_ID} + offset ${OFFSET}), exiting."
  exit 1
fi
 
echo "Task ${SLURM_ARRAY_TASK_ID} | Offset ${OFFSET} | Index ${ACTUAL_IDX} | Model: ${MODEL} | gene: ${GENE}"
 
Rscript run_gene_gemmamod_02.R "${GENE}"
 
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "ERROR: Rscript failed for gene ${GENE} with exit code ${EXIT_CODE}"
  exit $EXIT_CODE
fi
 
echo "Completed: ${GENE}"