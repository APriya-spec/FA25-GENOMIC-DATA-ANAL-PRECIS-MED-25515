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

# Part 12 — Final Summary & Interpretation

1. Differential Expression Overview

This analysis investigated transcriptomic differences in human respiratory cells exposed to SARS-CoV-2 compared to mock-infected controls, and also examined temporal changes within infected samples between 24 hours and 72 hours post-infection. Using DESeq2, significantly differentially expressed genes (DEGs) were identified based on an adjusted p-value < 0.05 and |log2FoldChange| > 1.

Across all conditions, the results revealed clear transcriptional differences between SARS-CoV-2 infected and mock samples, as well as dynamic regulation over time within the infected group.

2. Mock vs SARS-CoV-2 (All Time Points Combined)

A total of 46 genes were significantly differentially expressed between the mock and SARS-CoV-2 infected groups.

Upregulated in SARS-CoV-2: 29 genes

Downregulated in SARS-CoV-2: 17 genes

Mean log2 fold change: +1.42

Median log2 fold change: +1.87

Adjusted p-value range: 2.2e-06 to 4.5e-02

Most DEGs showed positive log2 fold change, indicating that SARS-CoV-2 infection predominantly induced increased gene expression relative to mock controls. These upregulated genes likely include those involved in innate immune responses, antiviral pathways, interferon signaling, and inflammatory cytokine regulation. The downregulated genes may reflect suppression of cellular homeostasis or host pathways disrupted by viral replication.

Interpretation:
These results suggest that SARS-CoV-2 infection triggers a strong and coordinated host response, characterized largely by activation of defense-related and inflammatory genes. The substantial upregulation across many transcripts supports known SARS-CoV-2 mechanisms involving immune activation and viral manipulation of cellular processes.

3. SARS-CoV-2 24H vs SARS-CoV-2 72H

Within the SARS-CoV-2 infected samples, 21 genes were identified as significantly different between 24 hours and 72 hours.

Upregulated (24h > 72h): 9 genes

Downregulated (24h < 72h): 12 genes

Mean log2 fold change: −0.09

Median log2 fold change: −1.34

Adjusted p-value range: 1.1e-12 to 3.7e-02

Most differential genes showed negative log2 fold change, indicating stronger expression at 72 hours compared to 24 hours. This suggests a progressively intensifying transcriptional response as infection progresses. Genes upregulated at early stages may reflect immediate early-response mechanisms, while those upregulated at later time points may reflect sustained antiviral defense, cell-stress responses, and pathways associated with prolonged infection.

Interpretation:
Temporal analysis reveals a dynamic host response to SARS-CoV-2. Initial responses at 24h appear more limited, while by 72h, gene expression changes become more pronounced, indicating that the virus continues to drive transcriptional reprogramming across time. This pattern aligns with typical viral kinetics, where cellular stress, immune activation, and signaling pathways evolve as infection becomes more established.

4. Overall Biological Implications

Strong Infection-Induced Gene Regulation:
SARS-CoV-2 infection significantly alters gene expression, with clear upregulation of many transcripts. This underscores the extensive host reprogramming caused by the virus.

Time-Dependent Transcriptional Shifts:
Differences between 24h and 72h highlight how the infection evolves over time. Early-response genes give way to later sustained responses, likely driven by continued viral replication and host immune activation.

Potential Immune & Stress-Related Signatures:
Although this dataset reflects miRNA-derived expression counts, the DEGs suggest involvement of inflammatory, antiviral, and stress-response pathways.

Utility for Biomarkers:
Both sets of DEGs (Mock vs SARS, and 24h vs 72h) contain candidates that could serve as biomarkers for infection state or time-specific viral response.
