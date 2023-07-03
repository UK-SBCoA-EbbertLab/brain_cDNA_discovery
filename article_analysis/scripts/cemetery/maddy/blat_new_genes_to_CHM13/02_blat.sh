#!/bin/bash
#SBATCH --time=40:15:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=pblat             # Job name
#SBATCH --ntasks=5                  # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=200G                  # Total memory requested
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err             # Error file for this job.
#SBATCH -o slurm-%j.out             # Output file for this job.
#SBATCH -A cca_mteb223_uksr         # Project allocation account name (REQUIRED)

singularity exec /home/mewa283/scripts/epigenetics_software_2023_02_22.sif pblat /project/mteb223_uksr/sequencing_resources/references/CHM13_plus_HG002_Y_chrom/GCA_009914755.4_CHM13_T2T_v2.0_genomic-RENAMED.fna $1 -out=blast9 $1.blast9
#blat /project/mteb223_uksr/sequencing_resources/references/CHM13_plus_HG002_Y_chrom/GCA_009914755.4_CHM13_T2T_v2.0_genomic-RENAMED.fna new_transcripts_for_new_genes.fa -out=blast9 output_new_transcripts_for_new_genes.blast9

