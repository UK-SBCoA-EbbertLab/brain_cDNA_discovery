#!/bin/bash


for fasta in /pscratch/mteb223_uksr/mlpa241/blat_new_genes_to_CHM13/new_genes_by_chr/*_new_transcripts_for_new_genes_exons.gtf.fa
do
	sbatch 02_blat.sh $fasta 

done
