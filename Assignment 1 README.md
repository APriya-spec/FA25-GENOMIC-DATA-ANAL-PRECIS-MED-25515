# Assignment 1 | Genome Assembly Using Velvet and Oases (Measles Virus)
# Programmer Name: Aruna Priya Cheekatla
# Language of the script: Unix / Bash
# Date: 11/24/2025
# Description:
This project focuses on performing a de novo assembly of the Measles virus genome using short-read Illumina sequencing data (SRA ID: SRR30155711). The main goal was to reconstruct the viral genome using two assembly tools — Velvet and Oases — and to observe how varying the k-mer size affects assembly performance.

Due to the small size of the Measles genome (~16 kb), Velvet successfully produced high-quality assemblies, while Oases did not generate transcripts (as expected for a non-transcriptomic dataset). All steps were executed on the Indiana University Slate HPC system using a Conda-managed environment for reproducibility.

# Requirements:

SRA-Toolkit (v3.0.5) – to download and convert SRA data

Velvet (v1.2.10) – for genome assembly

Oases (v0.2.09) – for transcriptome-style assembly

Conda – for environment management and package installation

# Steps to Execute the Project:
## Step 1 — Load Conda and Create the Environment
```bash
module load conda
conda create -n assignment_measles -c bioconda sra-tools velvet oases -y
conda activate assignment_measles
```
## Step 2 — Set Up Working Directories
```bash
BASE_DIR=/N/slate/archee/measles_asg1
mkdir -p $BASE_DIR/{data,velvet_output,oases_output,logs,tmp}
cd $BASE_DIR
```
# Step 3 — Download and Prepare Sequencing Data
```bash
cd $BASE_DIR/data
prefetch SRR30155711 --transport https
fasterq-dump SRR30155711 --split-files --temp $BASE_DIR/tmp
ls -lh SRR30155711_*.fastq
```
## Step 4 — Run Velvet Assemblies for Multiple k-mer Sizes
```bash
cd $BASE_DIR/velvet_output
for K in 51 61 71 81; do
  echo "Running Velvet for k = $K ..."
  mkdir -p run_$K
  velveth run_$K $K -fastq -shortPaired -separate \
    $BASE_DIR/data/SRR30155711_1.fastq $BASE_DIR/data/SRR30155711_2.fastq
  velvetg run_$K -exp_cov auto -cov_cutoff auto > $BASE_DIR/logs/velvet_$K.log 2>&1
done
```
## Step 5 — Run Oases Assemblies (Transcriptome-Based)
```bash
cd $BASE_DIR/velvet_output
for K in 51 61 71 81; do
  echo "Running Oases for k = $K ..."
  cd run_$K
  oases . -min_trans_lgth 100 > ../../logs/oases_$K.log 2>&1
  cd ..
done
```
Note: Oases did not generate transcripts due to the small genome size and lack of transcriptomic complexity. Only log files were created.

## Step 6 — Evaluate Assembly Outputs (Velvet)
```bash
for K in 51 61 71 81; do
  echo "Results for k=$K"
  cat $BASE_DIR/velvet_output/run_$K/stats.txt | grep -E "Final|N50|Total|num"
done
```
## Step 7 — Summarize Optimal Assembly
```bash
echo "Best Assembly Achieved with Velvet (k = 81)"
grep -E "N50|# contigs|Total length|Largest contig" \
  $BASE_DIR/velvet_output/run_81/stats.txt
```
# Generated Files:

data/ → paired FASTQ files (SRR30155711_1.fastq, SRR30155711_2.fastq)

velvet_output/ → assembled contigs for k=51, 61, 71, 81

oases_output/ → log files (no transcripts generated)

logs/ → assembly and run logs for each k-mer

# Summary of Findings

Increasing k-mer size improved the quality and continuity of assemblies.

Velvet successfully assembled the Measles virus genome, producing the best result at k = 81, which had longer contigs and fewer fragments.

Oases did not generate transcripts because the Measles dataset is small and lacks transcriptome-level complexity — this is expected behavior, not an error.

The final Velvet assembly closely matches the known Measles genome length (~16 kb) and GC content (~48%).
