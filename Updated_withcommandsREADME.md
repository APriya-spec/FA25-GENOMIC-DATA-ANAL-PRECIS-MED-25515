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

## Step 1 — Load Conda and Create the Environment
  
# Load the conda module on Slate HPC

module load conda

# Create a new conda environment for this assignment

conda create -n assignment_1_precision sra-tools velvet oases -y

# Activate the environment

conda activate assignment_1_precision

## Step 2 — Set Up Working Directories

# Define the base directory 

BASE_DIR=/N/slate/archee/ecoli_asg1

# Create folders for each stage of the project

mkdir -p $BASE_DIR/{data,velvet_output,oases_output,quast_results,logs,tmp}

# Navigate into the main directory

cd $BASE_DIR

## Step 3 — Download and Prepare Sequencing Data

# Move to the data folder

cd $BASE_DIR/data

# Download the E. coli short-read data from NCBI SRA

prefetch SRR21904868 --output-directory .

# Convert SRA file to paired FASTQ files (forward/reverse)

fasterq-dump SRR21904868 --split-files --temp $BASE_DIR/tmp

# List the FASTQ files to confirm

ls -lh SRR21904868_*.fastq

## Step 4 — Run Velvet Assemblies for Multiple k-mer Sizes

# Move to the velvet output folder

cd $BASE_DIR/velvet_output

# Run Velvet for k-mers 51, 61, 71, and 81

for K in 51 61 71 81; do
    echo "Running Velvet assembly for k = $K ..."
    mkdir -p run_$K
    velveth run_$K $K -fastq -shortPaired -separate \
        $BASE_DIR/data/SRR21904868_1.fastq $BASE_DIR/data/SRR21904868_2.fastq
    velvetg run_$K -exp_cov auto -cov_cutoff auto > $BASE_DIR/logs/velvet_$K.log 2>&1
done

## Step 5 — Run Oases Assemblies for the Same k-mer Sizes

# Move to the oases output folder

cd $BASE_DIR/oases_output

# Run Oases for k-mers 51, 61, 71, and 81

for K in 51 61 71 81; do
    echo "Running Oases assembly for k = $K ..."
    mkdir -p run_$K
    velveth run_$K $K -fastq -shortPaired -separate \
        $BASE_DIR/data/SRR21904868_1.fastq $BASE_DIR/data/SRR21904868_2.fastq
    velvetg run_$K -exp_cov auto -cov_cutoff auto
    oases run_$K > $BASE_DIR/logs/oases_$K.log 2>&1
done

### Step 6 — Evaluate Assemblies with QUAST

# Move back to base directory

cd $BASE_DIR

# Run QUAST for Velvet assemblies

quast velvet_output/run_*/contigs.fa -o quast_results/velvet_summary --min-contig 200

# Run QUAST for Oases assemblies

quast oases_output/run_*/transcripts.fa -o quast_results/oases_summary --min-contig 200

### Step 7 — View QUAST Reports

# View summary metrics for Velvet

less quast_results/velvet_summary/report.txt

# View summary metrics for Oases

less quast_results/oases_summary/report.txt

### Step 8 — Summarize Best Assembly

echo "Optimal Assembly: Velvet (k = 81)"

grep -E "N50|# contigs|Total length|Largest contig" \

 $BASE_DIR/quast_results/velvet_summary/report.txt



  ### Generated Files:

- data/ → contains paired FASTQ reads (SRR21904868_1.fastq, SRR21904868_2.fastq)

- velvet_output/ → contains contigs for each k-mer (51, 61, 71, 81)

- oases_output/ → contains transcripts for each k-mer (51, 61, 71, 81)

- quast_results/ → includes evaluation reports:

  Velvet QUAST Report

  Oases QUAST Report

- logs/ → contains assembly log files for each k-mer run

- assembly_pipeline.sh → script containing all commands used for this workflow

### Summary of Findings:

- Increasing k-mer size improved contiguity and assembly quality.

- Velvet (k = 81) generated the most accurate assembly (closest to 4.78 Mb), while smaller k-mers produced fragmented results.

- Oases (k = 81) produced longer transcripts but redundant contigs due to transcriptome-style assembly.

- Both assemblers maintained consistent GC content around 50%, matching E. coli reference values.


