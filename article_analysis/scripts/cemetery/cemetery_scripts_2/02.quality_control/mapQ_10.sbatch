#!/bin/bash
#SBATCH --time=5:00:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=flagstat_for_mapq10_%j    # Job name
#SBATCH -n 5 ## Ask for 31 cores
#SBATCH -N 1 ## Execute everything in only one node
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=normal         # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e flagstat_for_mapq10_%j.err             # Error file for this job.
#SBATCH -o flagstat_for_mapq10_%j.out             # Output file for this job.
#SBATCH -A cca_mteb223_uksr     # Project allocation account name (REQUIRED)

## Create flagstats for mapped reads with MAPQ > 10

for file_path in ../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_*filtered_sorted.bam
do
    sample_name="${file_path##*/}"

    singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif samtools flagstat \
        $file_path > "../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/${sample_name}.flagstat"

done
