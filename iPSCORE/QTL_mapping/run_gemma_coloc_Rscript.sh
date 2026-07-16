#!/bin/bash
#SBATCH --account=pro00117530
#SBATCH --job-name=gemma_coloc
#SBATCH --output=gemma_coloc.out
#SBATCH --error=gemma_coloc.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=48:00:00

# array-size needs to be specified within the sbatch command

module load miniforge
conda activate gemma_r422

Rscript GEMMA_res_coloc_Rscript.R
 
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "ERROR: Rscript failed for GEMMA results with exit code ${EXIT_CODE}"
  exit $EXIT_CODE
fi
 
echo "Completed GEMMA colocalization"