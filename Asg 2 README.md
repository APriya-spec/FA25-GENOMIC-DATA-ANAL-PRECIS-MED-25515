# Assignment 2 | miRNA Expression Analysis in SARS-CoV-2 Infected Cells
##  Programmer Name: Aruna Priya Cheekatla

##  Language of the scripts: Unix / Bash (alignment, QC, counting) and R (DESeq2 analysis)
##  Date: 11/16/2025

# Overview

This project performs transcriptomic analysis of human respiratory cells infected with SARS-CoV-2 compared to mock-treated cells at two different time points (24 hours and 72 hours). The aim is to identify differentially expressed microRNA genes and biological pathways altered during infection. The workflow includes downloading sequencing reads, quality control, trimming, alignment to the GRCh38 human genome, quantifying gene expression counts, and performing differential expression and GO enrichment analysis. The overall analysis provides insight into miRNA-mediated host responses to SARS-CoV-2 infection and temporal changes in gene regulation throughout the infection process.

# Requirements

## Software used:

- Cutadapt v4.9 – adapter & quality trimming

- HISAT2 v2.2.1 – genome indexing and alignment

- Samtools v1.17 – SAM/BAM manipulation

- Subread v2.0.3 (featureCounts) – read quantification

## R v4.4.2 with:

- DESeq2

- clusterProfiler

- org.Hs.eg.db

- ggplot2

- pheatmap

# PART 1 – Environment Setup

module load conda
conda create -n hisat2_env hisat2 fastqc trimgalore subread sra-tools samtools -y
conda activate hisat2_env

A dedicated environment ensures correct versions of HISAT2, FastQC, TrimGalore, and featureCounts are installed, making the workflow reproducible and isolated from system-wide changes.

# PART 2 – Directory Structure

mkdir -p /N/slate/archee/Janga_Asg2/{SRR,qc_raw,trimmed,reference,aligned,counts}
cd /N/slate/archee/Janga_Asg2

Organizing the project into separate folders for raw data, QC outputs, trimmed reads, reference files, alignments, and counts keeps the workflow clean and prevents file overwrites.

# PART 3 – Download Raw Sequencing Data

cd SRR
for id in SRR22269872 SRR22269873 SRR22269874 SRR22269875 SRR22269876 \
          SRR22269877 SRR22269878 SRR22269879 SRR22269880 SRR22269881 \
          SRR22269882 SRR22269883
do
    prefetch $id
    fasterq-dump $id -O /N/slate/archee/Janga_Asg2/SRR
done

All 12 SRA accessions for SARS-CoV-2 and mock samples are downloaded and converted to FASTQ format for downstream analysis.

# PART 4 – Quality Control with FastQC

cd /N/slate/archee/Janga_Asg2
fastqc SRR/*.fastq -o qc_raw

FastQC reports reveal adapter contamination, base quality issues, and sequence duplication levels, helping confirm the reads are suitable for trimming and alignment.

# PART 5 – TrimGalore Trimming (PHRED33)

cd trimmed
for fq in /N/slate/archee/Janga_Asg2/SRR/*.fastq
do
    trim_galore -q 30 --phred33 --length 16 --fastqc $fq -o /N/slate/archee/Janga_Asg2/trimmed
done

Trimming removes low-quality bases, adapters, and extremely short fragments, improving alignment efficiency and reducing mapping bias.

# PART 6 – Download GRCh38 Reference Genome (Ensembl)

cd /N/slate/archee/Janga_Asg2/reference

wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz

The human genome FASTA and GTF annotation provide the necessary reference for HISAT2 alignment and featureCounts gene quantification.

# PART 7 – Build HISAT2 Genome Index

hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa grch38_index

Indexing converts the FASTA reference into a searchable database, enabling fast and accurate read alignment using HISAT2.

# PART 8 – Alignment of Trimmed FASTQ Files

cd /N/slate/archee/Janga_Asg2/trimmed

for id in SRR22269872 SRR22269873 SRR22269874 SRR22269875 SRR22269876 \
          SRR22269877 SRR22269878 SRR22269879 SRR22269880 SRR22269881 \
          SRR22269882 SRR22269883
do
    hisat2 -p 8 \
        -x /N/slate/archee/Janga_Asg2/reference/grch38_index \
        -U ${id}_trimmed.fq.gz \
        -S /N/slate/archee/Janga_Asg2/aligned/${id}.sam \
        --rna-strandness R \
        --summary-file /N/slate/archee/Janga_Asg2/aligned/${id}_summary.txt
done

Reads map to the GRCh38 genome with strand-specific alignment. Summary reports show alignment rates between 47%–83%, indicating good read quality after trimming.

# PART 9 – Convert SAM to BAM

cd /N/slate/archee/Janga_Asg2/aligned

for id in SRR22269872 SRR22269873 SRR22269874 SRR22269875 SRR22269876 \
          SRR22269877 SRR22269878 SRR22269879 SRR22269880 SRR22269881 \
          SRR22269882 SRR22269883
do
    samtools view -S -b ${id}.sam > ${id}.bam
done

BAM format is compressed and required for featureCounts. Converting SAM to BAM reduces storage size and improves processing speed.

# PART 10 – Gene Quantification Using featureCounts

cd /N/slate/archee/Janga_Asg2

featureCounts \
  -a reference/Homo_sapiens.GRCh38.115.gtf \
  -o counts/gene_counts.txt \
  aligned/*.bam \
  -T 8 \
  -g gene_id

featureCounts assigns aligned reads to gene features. The resulting gene_counts.txt file becomes the input for DESeq2 differential expression analysis.

# PART 11 – Differential Expression Analysis (R Script)

The R code for this part is given in below link.

https://github.com/APriya-spec/FA25-GENOMIC-DATA-ANAL-PRECIS-MED-25515/blob/Asg-2/Janga%20asg2.R


