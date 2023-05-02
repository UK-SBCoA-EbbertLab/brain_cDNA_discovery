#!/bin/bash
#SBATCH --time=1:00:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=mafft_alignment_%j    # Job name
#SBATCH -n 1 ## Ask for 1 core
#SBATCH -N 1 ## Execute everything in only one node
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=normal         # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e mafft_alignment_%j.err             # Error file for this job.
#SBATCH -o mafft_alignment_%j.out             # Output file for this job.
#SBATCH -A cca_mteb223_uksr     # Project allocation account name (REQUIRED)


input_fasta=$1
output=$2

singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif mafft --auto --maxiterate 10000 $input_fasta > $output
