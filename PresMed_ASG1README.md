# **Assignment 1 | Genome Assembly Using Velvet and Oases**
### **Programmer Name: Aruna Priya Cheekatla**

Language of the script: Unix / Bash

Date: 10/14/2025

### **Description:**

This project focuses on the de novo assembly of the Escherichia coli genome using Illumina short-read sequencing data (SRA ID: SRR21904868). The main objective was to reconstruct the E. coli genome using two assembly tools — Velvet and Oases — and to study how varying k-mer sizes influence assembly quality. The quality of each assembly was assessed using QUAST, considering only contigs of at least 200 base pairs to eliminate small, low-quality fragments.

All steps were performed on Indiana University’s Slate HPC system, using a Conda environment to manage the installation of bioinformatics tools.

### **Requirements:**

1. SRA-Toolkit – for downloading and converting data from NCBI (version: 3.1.1)

2. Velvet – for genome assembly from short reads (version: 1.2.10)

3. Oases – transcriptome assembler based on Velvet (version: 0.2.09)

4. QUAST – for assembly quality evaluation (version: 5.2.0)

5. Conda – environment manager used to install and maintain software dependencies

### **Steps to Execute the Project:**

1. Logged into the Slate HPC environment and loaded the Conda module.

2. Created a Conda environment named assignment_1_precision with sra-tools, velvet, and oases.

3. Created a project directory named ecoli_asg1 with subfolders for data, outputs, logs, and results.

4. Downloaded the Illumina short-read data (SRR21904868) from NCBI and converted it to paired FASTQ files.

5. Performed Velvet assemblies using k-mer values of 51, 61, 71, and 81.

6. Ran Oases assemblies with the same k-mer values to compare transcriptome-based performance.

7. Evaluated both assemblies using QUAST with a minimum contig length of 200 bp to remove small and unreliable fragments.

8. Compared assembly metrics such as total length, N50, GC%, and number of contigs to determine the best performing k-mer size and assembler.

### **Generated Files:**

1. Data folder: Contains paired-end FASTQ reads (SRR21904868_1.fastq and SRR21904868_2.fastq).

2. Velvet Output folder: Includes assembled contigs and statistics for each k-mer tested (51, 61, 71, 81).

3. Oases Output folder: Contains transcript assemblies and related results for each k-mer run.

4. QUAST Results: Contains detailed evaluation reports for Velvet and Oases assemblies.
   
https://github.com/APriya-spec/FA25-GENOMIC-DATA-ANAL-PRECIS-MED-25515/blob/Asg-1/oases_summary.zip

https://github.com/APriya-spec/FA25-GENOMIC-DATA-ANAL-PRECIS-MED-25515/blob/Asg-1/velvet_summary.zip

5. Logs folder: Includes run logs for both Velvet and Oases showing execution details and assembly progress.

### **Summary of Findings:**

Assemblies produced with smaller k-mers (51, 61) were highly fragmented and had lower N50 values, while larger k-mers (71, 81) resulted in more contiguous assemblies with longer contigs. The Velvet assembly at k = 81 yielded the best genome reconstruction (~4.78 Mb total length, N50 = 112 kb), closely matching the known E. coli genome size. The Oases assembly at k = 81 also showed good results but contained redundant sequences due to overlapping transcripts.

Using a 200 bp contig threshold improved assembly reliability and provided a clearer view of quality across k-mer values. Overall, Velvet at k = 81 produced the most accurate and complete E. coli genome assembly.
