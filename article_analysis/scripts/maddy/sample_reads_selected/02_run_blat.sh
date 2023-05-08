#!/bin/bash


for fasta in /pscratch/mteb223_uksr/brain_cDNA_discovery/article_analysis/data/maddy/sample_reads_selected/*_BambuGene290099.fa 
do
	sbatch /pscratch/mteb223_uksr/mlpa241/blat_new_genes_to_CHM13/02_blat.sh $fasta 

done
