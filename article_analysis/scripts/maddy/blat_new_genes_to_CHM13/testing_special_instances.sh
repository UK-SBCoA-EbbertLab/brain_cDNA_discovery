#!/bin/bash
#SBATCH --time=40:15:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=blat             # Job name
#SBATCH --ntasks=1                  # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=100G                  # Total memory requested
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err             # Error file for this job.
#SBATCH -o slurm-%j.out             # Output file for this job.
#SBATCH -A cca_mteb223_uksr         # Project allocation account name (REQUIRED)

blat /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa test_new_poster_child.fa -out=blast9 output_test_back_to_hg38.blast9

