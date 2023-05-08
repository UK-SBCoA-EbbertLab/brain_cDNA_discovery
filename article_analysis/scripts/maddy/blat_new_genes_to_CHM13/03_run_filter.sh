#!/bin/bash


for blast in /pscratch/mteb223_uksr/mlpa241/blat_new_genes_to_CHM13/new_genes_by_chr/*_new_transcripts_for_new_genes.gtf.fa.blast9
do
	singularity exec /project/mteb223_uksr/singularity_files/transcriptome_long_reads_GTEx_2022_12_14.sif Rscript 03_filter_blast_results.R $blast

done
