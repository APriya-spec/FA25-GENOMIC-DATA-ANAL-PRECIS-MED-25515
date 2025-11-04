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
module load conda

- Create and activate the environment:
conda create -n assignment_1_precision sra-tools velvet oases
conda activate assignment_1_precision

2. Directory Structure

- Created a project folder named ecoli_asg1 with subfolders:
data/, velvet_output/, oases_output/, quast_results/, and logs/.
