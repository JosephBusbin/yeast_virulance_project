#!/bin/bash
#SBATCH --job-name=armaan_one_line
#SBATCH --nodes=1 # how many physical machines in the cluster
#SBATCH --ntasks=1 # how many separate 'tasks' (stick to 1)
#SBATCH --cpus-per-task=6 # how many cores (bora max is 20)
#SBATCH --time=120 # d-hh:mm:ss or just No. of minutes
#SBATCH -o oneline_%j.out #STDOUT to file (%j is jobID)
#SBATCH -e oneline_%j.err #STDERR to file (%j is jobID)

set -ueo pipefail

# ==============================================================================
# Pipeline Variables
# Usage: ./run_pipeline.sh <SAMPLE_NAME> <OUTPUT_DIR>
# ==============================================================================
SAMPLE="${1:-sample1}"
OUT_DIR="${2:-output}"
REF_FASTA="reference.fasta"
FASTQ_R1="${SAMPLE}.fastq.gz"
THREADS=8

# Read Group metadata
RG_ID="${SAMPLE}"
RG_SM="${SAMPLE}"
RG_LB="Murphy_Lab"
RG_PL="illumina"
RG_PU="unit1"

# Directory Structure
RAW_VCF_DIR="${OUT_DIR}/raw_vcf"
FILTERED_VCF_DIR="${OUT_DIR}/filtered_vcf"
mkdir -p "${RAW_VCF_DIR}" "${FILTERED_VCF_DIR}"

# Intermediate BAM & Metric File Paths
BAM_MEM="${OUT_DIR}/${SAMPLE}.mem.bam"
BAM_SORTED="${OUT_DIR}/${SAMPLE}.sorted.mem.bam"
BAM_RG="${OUT_DIR}/${SAMPLE}.sorted.rg.bam"
BAM_DEDUP="${OUT_DIR}/${SAMPLE}.markdup.rg.bam"
DUP_METRICS="${OUT_DIR}/${SAMPLE}.dup_metrics.txt"

# Output VCF File Paths
RAW_VCF="${RAW_VCF_DIR}/${SAMPLE}.gatk.vcf"
FILTERED_VCF="${FILTERED_VCF_DIR}/${SAMPLE}.Filtered.vcf"
SUBSET_VCF="${FILTERED_VCF_DIR}/${SAMPLE}.FilteredSubset.vcf"

# Derived reference dictionary name required by GATK
REF_DICT="${REF_FASTA%.*}.dict"

# ==============================================================================
# Environment & Module Loading
# ==============================================================================
module load gatk/4.6.1.0
module load samtools/gcc-11.4.1/1.21
module load bwa/gcc-11.4.1/0.7.18

# ==============================================================================
# Step 0: Reference Genome Preparation (Executed only if indices are missing)
# ==============================================================================
if [ ! -f "${REF_DICT}" ]; then
    gatk CreateSequenceDictionary -R "${REF_FASTA}" -O "${REF_DICT}"
fi

if [ ! -f "${REF_FASTA}.fai" ]; then
    samtools faidx "${REF_FASTA}"
fi

if [ ! -f "${REF_FASTA}.bwt" ]; then
    bwa index "${REF_FASTA}"
fi

# ==============================================================================
# Step 1: Align Single-End Reads
# ==============================================================================
bwa mem -t "${THREADS}" "${REF_FASTA}" "${FASTQ_R1}" | samtools view -bS - > "${BAM_MEM}"

# ==============================================================================
# Step 2: Sort Alignment by Genomic Coordinates
# ==============================================================================
samtools sort -@ "${THREADS}" "${BAM_MEM}" -o "${BAM_SORTED}"

# ==============================================================================
# Step 3: Add Read Groups (GATK)
# ==============================================================================
gatk AddOrReplaceReadGroups \
    -I "${BAM_SORTED}" \
    -O "${BAM_RG}" \
    -RGPL "${RG_PL}" \
    -RGID "${RG_ID}" \
    -RGPU "${RG_PU}" \
    -RGLB "${RG_LB}" \
    -RGSM "${RG_SM}"

# ==============================================================================
# Step 4: Mark Duplicates (GATK)
# ==============================================================================
gatk MarkDuplicates \
    -I "${BAM_RG}" \
    -O "${BAM_DEDUP}" \
    -M "${DUP_METRICS}" \
    --REMOVE_DUPLICATES true

# ==============================================================================
# Step 5: Index Deduplicated BAM (GATK)
# ==============================================================================
gatk BuildBamIndex -I "${BAM_DEDUP}"

# ==============================================================================
# Step 6: Variant Calling (HaplotypeCaller) -> Saves to raw_vcf/
# ==============================================================================
gatk HaplotypeCaller \
    -R "${REF_FASTA}" \
    -I "${BAM_DEDUP}" \
    -O "${RAW_VCF}"

# ==============================================================================
# Cleanup: Purge Non-VCF Intermediate Files
# ==============================================================================
rm -f "${BAM_MEM}" \
      "${BAM_SORTED}" \
      "${BAM_RG}" \
      "${BAM_DEDUP}" \
      "${BAM_DEDUP}.bai" \
      "${OUT_DIR}/${SAMPLE}.markdup.rg.bai" \
      "${DUP_METRICS}"

# ==============================================================================
# Step 7: Filter Variants with Hard Filters -> Saves to filtered_vcf/
# ==============================================================================
gatk VariantFiltration \
    -R "${REF_FASTA}" \
    -V "${RAW_VCF}" \
    -O "${FILTERED_VCF}" \
    --filter-name "snpsFilter1" \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 2.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

# ==============================================================================
# Step 8: Subset Filtered Variants (PASS only) -> Saves to filtered_vcf/
# ==============================================================================
gatk SelectVariants \
    -R "${REF_FASTA}" \
    -V "${FILTERED_VCF}" \
    --exclude-filtered \
    -O "${SUBSET_VCF}"

# ==============================================================================
# Step 9: Validate Output VCF
# ==============================================================================
gatk ValidateVariants \
    -R "${REF_FASTA}" \
    -V "${SUBSET_VCF}"