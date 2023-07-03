#!/bin/bash
#SBATCH --time=01:15:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=minimap_index    # Job name
#SBATCH --ntasks=8                  # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --nodes 1
#SBATCH --mem=70G                  # Total memory requested
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err             # Error file for this job.
#SBATCH -o slurm-%j.out             # Output file for this job.
#SBATCH -A cca_mteb223_uksr         # Project allocation account name (REQUIRED)

singularity exec /pscratch/mteb223_uksr/brain_cDNA_discovery/singularity_containers/nanopore.sif minimap2 -t 8 -d chm13_ref.mmi /project/mteb223_uksr/sequencing_resources/references/CHM13_plus_HG002_Y_chrom/GCA_009914755.4_CHM13_T2T_v2.0_genomic-RENAMED.fna

