### Assignment 1 | Genome Assembly Using Velvet and Oases
### Programmer Name: Aruna Priya Cheekatla

### Language of the script: Unix / Bash
### Date: 10/14/2025
### Description:

This project involves de novo assembly of the Escherichia coli genome using short-read Illumina sequencing data (SRA ID: SRR21904868).
The main objective was to reconstruct the genome using two assembly tools — Velvet and Oases — and evaluate how varying k-mer sizes influence assembly quality.
The quality of each assembly was assessed using QUAST, considering only contigs longer than 200 base pairs to remove low-quality fragments.

All analyses were conducted on the Indiana University Slate HPC system using a Conda-managed environment.

### Requirements:

SRA-Toolkit – to download and convert raw sequencing data (v3.1.1)

Velvet – for de novo genome assembly (v1.2.10)

Oases – for transcriptome assembly (v0.2.09)

QUAST – for assessing assembly quality (v5.2.0)

Conda – environment manager for reproducibility

### Steps to Execute the Project:

1. Environment Setup

- Load the Conda module on Slate HPC:
- 
module load conda

- Create and activate the environment:
  
conda create -n assignment_1_precision sra-tools velvet oases

conda activate assignment_1_precision

2. Directory Structure

- Created a project folder named ecoli_asg1 with subfolders:
  
data/, velvet_output/, oases_output/, quast_results/, and logs/.

3. Data Download and Preparation

- Downloaded sequencing reads from NCBI SRA using:
  
prefetch SRR21904868

- Converted .sra to paired-end FASTQ files using:
  
fasterq-dump SRR21904868 --split-files

4. Genome Assembly with Velvet

- Performed assemblies for multiple k-mer sizes: 51, 61, 71, and 81.

- Example command:

velveth run_81 81 -fastq -shortPaired -separate SRR21904868_1.fastq SRR21904868_2.fastq

velvetg run_81 -exp_cov auto -cov_cutoff auto

5. Assembly with Oases

- Repeated the assembly process using Oases for the same k-mer values.

- Example command:

velveth run_81 81 -fastq -shortPaired -separate SRR21904868_1.fastq SRR21904868_2.fastq

velvetg run_81 -exp_cov auto -cov_cutoff auto

oases run_81

6. Assembly Quality Evaluation (QUAST)

- Compared all assemblies using QUAST with a minimum contig size of 200 bp:

quast velvet_output/run_*/contigs.fa -o quast_results/velvet_summary --min-contig 200

quast oases_output/run_*/transcripts.fa -o quast_results/oases_summary --min-contig 200

7. Result Comparison and Optimization

- Evaluated the QUAST reports to identify the optimal k-mer for each tool.

- Determined Velvet (k=81) produced the most contiguous assembly, while Oases performed better for transcript reconstruction but included redundant regions.

  ### Generated Files:

data/ → contains paired FASTQ reads (SRR21904868_1.fastq, SRR21904868_2.fastq)

velvet_output/ → contains contigs for each k-mer (51, 61, 71, 81)

oases_output/ → contains transcripts for each k-mer (51, 61, 71, 81)

quast_results/ → includes evaluation reports:

Velvet QUAST Report

Oases QUAST Report

logs/ → contains assembly log files for each k-mer run

assembly_pipeline.sh → script containing all commands used for this workflow

### Summary of Findings:

Increasing k-mer size improved contiguity and assembly quality.

Velvet (k = 81) generated the most accurate assembly (closest to 4.78 Mb), while smaller k-mers produced fragmented results.

Oases (k = 81) produced longer transcripts but redundant contigs due to transcriptome-style assembly.

Both assemblers maintained consistent GC content around 50%, matching E. coli reference values.


