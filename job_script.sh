#!/bin/bash
#SBATCH --job-name=armaan_one_line
#SBATCH --nodes=1 # how many physical machines in the cluster
#SBATCH --ntasks=1 # how many separate 'tasks' (stick to 1)
#SBATCH --cpus-per-task=6 # how many cores (bora max is 20)
#SBATCH --time=120 # d-hh:mm:ss or just No. of minutes
#SBATCH -o oneline_%j.out #STDOUT to file (%j is jobID)
#SBATCH -e oneline_%j.err #STDERR to file (%j is jobID)

set -ueo pipefail

#job_script_for_one_lines

# file path ./data/raw/B02-A8a2-2338_S5_L002_R1_001.fastq
# reference ./data/reference/Anc8.genome.fasta

# ==============================================================================
# Pipeline Variables (Override sample name via argument 1: ./run_pipeline.sh <SAMPLE>)
# ==============================================================================
SAMPLE="B02-A8a2-2338_S5_L002_R1_001"
#"${1:-sample1}"
REF_FASTA="./data/reference/Anc8.genome.fasta"
#"reference.fasta"
FASTQ_R1="./data/raw/B02-A8a2-2338_S5_L002_R1_001.fastq"
#"${SAMPLE}.fastq.gz" #or "${SAMPLE}.fastq"
THREADS=6 #Depends on batch script

# Read Group metadata
RG_ID="${SAMPLE}"
RG_SM="${SAMPLE}"
RG_LB="Murphy_Lab"
RG_PL="illumina"
RG_PU="unit1"

# Derived reference dictionary name required by GATK
REF_DICT="${REF_FASTA%.*}.dict"

# ==============================================================================
# Environment & Module Loading
# ==============================================================================
module load gatk/4.6.1.0
module load samtools/gcc-11.4.1/1.21
module load bwa/gcc-11.4.1/0.7.18

# Initialize and activate Conda environment safely in non-interactive shell
eval "$(conda shell.bash hook)"
conda activate yeast_env

# ==============================================================================
# Step 0: Reference Genome Preparation (Executed only if indices are missing)
# ==============================================================================
if [ ! -f "${REF_DICT}" ]; then
    picard CreateSequenceDictionary R="${REF_FASTA}" O="${REF_DICT}"
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
#bwa mem -t "${THREADS}" "${REF_FASTA}" "${FASTQ_R1}" | samtools view -bS - > "${SAMPLE}.mem.bam"

# ==============================================================================
# Step 2: Sort Alignment by Genomic Coordinates
# ==============================================================================
#samtools sort -@ "${THREADS}" "${SAMPLE}.mem.bam" -o "${SAMPLE}.sorted.mem.bam"

# ==============================================================================
# Step 3: Add Read Groups (Picard)
# ==============================================================================
#picard AddOrReplaceReadGroups \
    #I="${SAMPLE}.sorted.mem.bam" \
    #O="${SAMPLE}.sorted.rg.bam" \
    #RGPL="${RG_PL}" \
    #RGID="${RG_ID}" \
    #RGPU="${RG_PU}" \
    #RGLB="${RG_LB}" \
    #RGSM="${RG_SM}"

# ==============================================================================
# Step 4: Mark Duplicates
# ==============================================================================
#picard MarkDuplicates \
    #I="${SAMPLE}.sorted.rg.bam" \
    #O="${SAMPLE}.markdup.rg.bam" \
    #M="${SAMPLE}.dup_metrics.txt" \
    #REMOVE_DUPLICATES=true


# ==============================================================================
# Step 5: Index Deduplicated BAM
# ==============================================================================
#picard BuildBamIndex I="${SAMPLE}.markdup.rg.bam"

# ==============================================================================
# Step 6: Variant Calling (HaplotypeCaller)
# ==============================================================================
gatk HaplotypeCaller \
    -R "${REF_FASTA}" \
    -I "${SAMPLE}.markdup.rg.bam" \
    -O "${SAMPLE}.gatk.vcf"
