# Assignment 2 | miRNA Expression Analysis in SARS-CoV-2 Infected Cells
# Programmer Name: Aruna Priya Cheekatla

# Language of the scripts: Unix / Bash (alignment, QC, counting) and R (DESeq2 analysis)
# Date: 11/16/2025

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
