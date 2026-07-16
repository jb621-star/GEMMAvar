#!/bin/bash
#SBATCH --account=pro00117530
#SBATCH --job-name=gemmamod_coloc_array
#SBATCH --output=logs/gemmamod_coloc_%A_%a.out
#SBATCH --error=logs/gemmamod_coloc_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=48:00:00
#SBATCH --array=1-12    # adjust upper bound to n_chunks

mkdir -p logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gemma_r422   # or whichever env has R + coloc

CHUNK=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

Rscript run_gemmamod_coloc_chunk.R \
    --chunk_file "pair_chunks/pairs_chunk_${CHUNK}.csv" \
    --outfile    "coloc_results/gemmamod_coloc_results_${CHUNK}.csv"