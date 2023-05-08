#!/bin/bash

#gtf=$1

#awk -F '\t' '$3 =="transcript" && $9 ~ /BambuGene/ {print > $1"_new_transcripts_for_new_genes.gtf"}' '/pscratch/mteb223_uksr/brain_cDNA_discovery/article_analysis/data/bernardo/processed/99.other/create_annotations/annotations_and_quant_for_mark_and_maddy/unfiltered_novel_genes_and_transcripts/new_rna_UNFILTERED_cpm_greater_than_one.gtf'
awk -F '\t' '$3 =="exon" && $9 ~ /BambuGene/ {print > "new_genes_by_chr/"$1"_new_transcripts_for_new_genes_exons.gtf"}' '/pscratch/mteb223_uksr/brain_cDNA_discovery/article_analysis/data/bernardo/processed/99.other/create_annotations/annotations_and_quant_for_mark_and_maddy/unfiltered_novel_genes_and_transcripts/new_rna_UNFILTERED_cpm_greater_than_one.gtf'


# filtered by median CPM > 1
#awk -F '\t' '$3 =="exon" && $9 ~ /BambuGene/ {print > "new_genes_by_chr/"$1"_new_transcripts_for_new_genes_filtered_MEDIAN_CPM.gtf"}' '/pscratch/mteb223_uksr/brain_cDNA_discovery/article_analysis/data/bernardo/processed/99.other/create_annotations/annotations_and_quant_for_mark_and_maddy/filtered_novel_genes_and_transcripts/new_rna_MEDIAN_cpm_greater_than_one.gtf'
#awk -F '\t' '$3 =="transcript" && $9 ~ /BambuGene/' $1 > new_transcripts_for_new_genes.gtf

singularity exec /project/mteb223_uksr/singularity_files/transcriptome_long_reads_GTEx_2022_12_14.sif Rscript 01a_create_key.R 

#for file in new_genes_by_chr/*_new_transcripts_for_new_genes.gtf
for file in new_genes_by_chr/*_new_transcripts_for_new_genes_exons.gtf
do
	singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed $file -fo $file.fa
done
