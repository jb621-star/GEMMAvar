#!/bin/bash
#SBATCH --job-name=prep_1kg
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --array=1-22:1
#SBATCH --cpus-per-task=8
#SBATCH --export=ALL
#SBATCH --partition=scavenger
#SBATCH --mem=8G

module load Plink/2
module load bcftools

pca_dir=/work/jb621/GEMMAmod_iPSCORE/prepare_wgs/pca
chr=$SLURM_ARRAY_TASK_ID

pull_vcf=https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget $pull_vcf 

ref_vcf=ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

# 1. Separate multi-allelic snps
cmd="bcftools norm --threads 8 -m - -Oz -o ${pca_dir}/chr${chr}.vcf.gz $ref_vcf"
echo $cmd; eval $cmd

# 2. Extract SNPs only, rename variants, output as PLINK bed files
cmd="plink2 --snps-only --vcf ${pca_dir}/chr${chr}.vcf.gz --make-bed --memory 3000 --threads 8 --out ${pca_dir}/chr${chr}.rename --set-all-var-ids @:#:\\\$1:\\\$2"
echo $cmd; eval $cmd

# 3. Write filename to merge list (with file locking to avoid race conditions)
cmd="echo ${pca_dir}/chr${chr}.rename"
echo $cmd
flock ${pca_dir}/merge_list.txt.lock -c "echo ${pca_dir}/chr${chr}.rename >> ${pca_dir}/merge_list.txt"
