#!/bin/bash
#SBATCH --time=1-00:00:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=extract_mito_bams_%j    # Job name
#SBATCH --mem=100G                  ## Ask for 512G of RAM
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e extract_mito_bams_%j.err             # Error file for this job.
#SBATCH -o extract_mito_bams_%j.out             # Output file for this job.
#SBATCH -A coa_jbmi230_uksr     # Project allocation account name (REQUIRED)



singularity exec ../../../../../singularity_containers/bernardo_article_analysis.sif ./extract_new_reads_to_bam.sh
