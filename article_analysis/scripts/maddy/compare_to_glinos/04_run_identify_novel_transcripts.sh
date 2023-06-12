#!/bin/bash
#SBATCH --time=24:15:00                                         # Time limit for the job (REQUIRED).
#SBATCH --job-name=GTF_comparisons      # Job name
#SBATCH --ntasks=1                                              # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=100G                                                # Total memory requested
#SBATCH --partition=normal                                      # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err                                         # Error file for this job.
#SBATCH -o slurm-%j.out                                         # Output file for this job.
#SBATCH -A coa_mteb223_uksr                                     # Project allocation account name (REQUIRED)


#their_data="Homo_sapiens.GRCh38.105.simplified.gtf"
#output_file="our_data_vs_Ensembl_condensed_100bp_threshold.tsv"

their_data="../../../data/maddy/compare_to_glinos/GTEx_GRCh38_version_107_mapq_10_track_reads/bambu_discovery/extended_annotations.gtf"
#their_data="Homo_sapiens.GRCh38.105.simplified.gtf"
#output_file_prefix="Glinos_vs_Ensembl_condensed_threshold"
#bp_threshold=100
novel_column_alt_value="GLINOS"

our_data=""
output_file="Ours_vs_Glinos_condensed_threshold.all.tsv"


#our_data="extended_annotation_grch38_only_novel.sorted.gtf"
#their_data="flair_filter_transcripts_no_chr.sorted.gtf"
#output_file="our_novel_transcripts_condensed_100bp_threshold.tsv"
#bp_threshold=100
#novel_column_alt_value="GLINOS_OVERLAP"



singularity exec /project/mteb223_uksr/singularity_files/transcriptome_long_reads_GTEx.sif Rscript 04_identify_novel_transcripts.R $our_data $their_data $output_file $novel_column_alt_value
