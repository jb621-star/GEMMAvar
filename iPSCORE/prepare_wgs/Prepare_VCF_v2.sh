#!/bin/bash
#SBATCH --account=pro00117530
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=12-00:00:00
#SBATCH --job-name=liftOver_hg38
#SBATCH --output=liftOver_hg38.out
#SBATCH --error=liftOver_hg38.err

module load miniforge
conda activate vcf_tools

# ── Paths ────────────────────────────────────────────────────────────────────
indir=/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/phg000904.v1.NHLBI_iPSCORE_WGS.genotype-calls-vcf.c1
outdir=/data/irb/biostatisticsbioinformatics/pro00117530/GEMMAmod_iPSCORE/prepare_wgs
#refdir=/data/irb/biostatisticsbioinformatics/pro00117530/ref

vcf_hg19=${indir}/IPSCORE_WGS_c1.chr.vcf.gz
chain=${indir}/hg19ToHg38.over.chain.gz
ref_hg38=${indir}/hg38.fa
log=${outdir}/liftOver_hg38.log

mkdir -p ${outdir}

# ── Step 1: Normalize hg19 VCF before liftover ───────────────────────────────
# Split multiallelic sites and left-align indels BEFORE lifting.
# This prevents CrossMap from producing malformed multiallelic records.
#echo "[$(date)] Step 1: Normalizing hg19 VCF" | tee -a $log

#cmd="bcftools norm --threads 8 \
#    -m - \
#    -f ${indir}/hg19.fa \
#    -Oz -o ${outdir}/IPSCORE_WGS_c1.norm.vcf.gz \
#    ${vcf_hg19}"
#echo $cmd | tee -a $log; eval $cmd

#cmd="tabix -p vcf ${outdir}/IPSCORE_WGS_c1.norm.vcf.gz"
#echo $cmd | tee -a $log; eval $cmd

# ── Step 2: Run CrossMap ──────────────────────────────────────────────────────
# CrossMap outputs:
#   *.hg38.vcf.gz        → successfully lifted variants
#   *.hg38.vcf.gz.unmap  → variants that failed to lift (inspect these!)
#echo "[$(date)] Step 2: Running CrossMap" | tee -a $log

#cmd="CrossMap vcf \
#    --chromid s \
#    ${chain} \
#    ${outdir}/IPSCORE_WGS_c1.norm.vcf.gz \
#    ${ref_hg38} \
#    ${outdir}/IPSCORE_WGS_c1.hg38.vcf"
#echo $cmd | tee -a $log; eval $cmd

# ── Step 3: Sort the lifted VCF ──────────────────────────────────────────────
# CrossMap does NOT guarantee coordinate order after lifting.
# An unsorted VCF will break all downstream tools.
#echo "[$(date)] Step 3: Sorting lifted VCF" | tee -a $log

#cmd="bcftools sort \
#    --temp-dir ${outdir}/tmp \
#    -Oz -o ${outdir}/IPSCORE_WGS_c1.hg38.sorted.vcf.gz \
#    ${outdir}/IPSCORE_WGS_c1.hg38.vcf"
#echo $cmd | tee -a $log; eval $cmd

#cmd="tabix -p vcf ${outdir}/IPSCORE_WGS_c1.hg38.sorted.vcf.gz"
#echo $cmd | tee -a $log; eval $cmd

# ── Step 4: Re-normalize against hg38 reference ──────────────────────────────
echo "[$(date)] Step 4: Re-normalizing against hg38" | tee -a $log

canonical="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

cmd="bcftools norm --threads 8 \
    -m - \
    -f ${ref_hg38} \
    --check-ref ws \
    --regions ${canonical} \
    -Oz -o ${outdir}/IPSCORE_WGS_c1.hg38.norm.vcf.gz \
    ${outdir}/IPSCORE_WGS_c1.hg38.sorted.chrfix.vcf.gz"
echo $cmd | tee -a $log; eval $cmd

cmd="tabix -p vcf ${outdir}/IPSCORE_WGS_c1.hg38.norm.vcf.gz"
echo $cmd | tee -a $log; eval $cmd

# ── Step 5: Remove problem variants ──────────────────────────────────────────
# Exclude: symbolic alleles (<DEL>, <NON_REF>, etc.), spanning deletions (*),
# duplicate positions, and non-SNP/indel records.
echo "[$(date)] Step 5: Cleaning problematic variants" | tee -a $log

cmd="bcftools view --threads 8 \
    -e 'ALT~\"<\" || ALT=\"*\" || ALT=\".\"' \
    ${outdir}/IPSCORE_WGS_c1.hg38.norm.vcf.gz \
    | bcftools view --threads 8 -v snps,indels \
    | bcftools norm --threads 8 -d exact \
    -Oz -o ${outdir}/IPSCORE_WGS_c1.hg38.norm.clean.vcf.gz"
echo $cmd | tee -a $log; eval $cmd

cmd="tabix -p vcf ${outdir}/IPSCORE_WGS_c1.hg38.norm.clean.vcf.gz"
echo $cmd | tee -a $log; eval $cmd

# ── Step 6: Filter with VCFtools ─────────────────────────────────────────────
echo "[$(date)] Step 6: Filtering with VCFtools" | tee -a $log

cmd="vcftools \
    --gzvcf ${outdir}/IPSCORE_WGS_c1.hg38.norm.clean.vcf.gz \
    --maf 0.05 \
    --max-alleles 2 \
    --hwe 1e-06 \
    --max-missing 0.99 \
    --remove-filtered-all \
    --recode --recode-INFO-all \
    --out ${outdir}/IPSCORE_WGS_c1.hg38.norm.clean.filtered"
echo $cmd | tee -a $log; eval $cmd

# ── Step 7: Bgzip + index the filtered VCF ───────────────────────────────────
echo "[$(date)] Step 7: Compressing and indexing" | tee -a $log

cmd="bgzip -@ 8 -c \
    ${outdir}/IPSCORE_WGS_c1.hg38.norm.clean.filtered.recode.vcf \
    > ${outdir}/IPSCORE_WGS_c1.hg38.norm.clean.filtered.vcf.gz"
echo $cmd | tee -a $log; eval $cmd

cmd="tabix -p vcf ${outdir}/IPSCORE_WGS_c1.hg38.norm.clean.filtered.vcf.gz"
echo $cmd | tee -a $log; eval $cmd

# ── Step 8: Rename variant IDs to CHROM_POS_REF_ALT ─────────────────────────
echo "[$(date)] Step 8: Renaming variant IDs" | tee -a $log

cmd="bcftools annotate --threads 8 \
    --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    ${outdir}/IPSCORE_WGS_c1.hg38.norm.clean.filtered.vcf.gz \
    -Oz -o ${outdir}/IPSCORE_WGS_c1.hg38.final.vcf.gz"
echo $cmd | tee -a $log; eval $cmd

cmd="tabix -p vcf ${outdir}/IPSCORE_WGS_c1.hg38.final.vcf.gz"
echo $cmd | tee -a $log; eval $cmd

# ── QC summary ───────────────────────────────────────────────────────────────
echo "[$(date)] Done. Variant counts:" | tee -a $log
echo -n "  hg19 input:    " | tee -a $log
bcftools stats ${vcf_hg19} | grep "^SN" | grep "number of records" | tee -a $log
echo -n "  post-liftover: " | tee -a $log
bcftools stats ${outdir}/IPSCORE_WGS_c1.hg38.norm.clean.vcf.gz | grep "^SN" | grep "number of records" | tee -a $log
echo -n "  post-filter:   " | tee -a $log
bcftools stats ${outdir}/IPSCORE_WGS_c1.hg38.final.vcf.gz | grep "^SN" | grep "number of records" | tee -a $log
echo -n "  unlifted:      " | tee -a $log
wc -l < ${outdir}/IPSCORE_WGS_c1.hg38.vcf.gz.unmap | tee -a $log