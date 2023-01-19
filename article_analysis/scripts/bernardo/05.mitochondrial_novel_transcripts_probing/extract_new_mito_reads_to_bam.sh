#!/bin/bash

## Only keep the first column that contains read ids and skip the header
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/PAM54401_1271_nanopore_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "uky_1271_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/PAM54902_1291_nanopore_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "uky_1291_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/PAM54788_1304_nanopore_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "uky_1304_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/PAM54335_356_nanopore_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "uky_356_mito_read_ids.txt"}'

## Extract the read ids mapping to new mitochondrial transcripts from the BAM files
samtools view ../../../data/bernardo/raw/uky_aged_stringent/minimap2/PAM54401_1271_nanopore_mapped_filtered_sorted.bam | fgrep -w -f uky_1271_mito_read_ids.txt > new_mito_reads_uky_1271.sam
samtools view ../../../data/bernardo/raw/uky_aged_stringent/minimap2/PAM54902_1291_nanopore_mapped_filtered_sorted.bam | fgrep -w -f uky_1291_mito_read_ids.txt > new_mito_reads_uky_1291.sam
samtools view ../../../data/bernardo/raw/uky_aged_stringent/minimap2/PAM54788_1304_nanopore_mapped_filtered_sorted.bam | fgrep -w -f uky_1304_mito_read_ids.txt > new_mito_reads_uky_1304.sam
samtools view ../../../data/bernardo/raw/uky_aged_stringent/minimap2/PAM54335_356_nanopore_mapped_filtered_sorted.bam | fgrep -w -f uky_356_mito_read_ids.txt > new_mito_reads_uky_356.sam


## Get header from bam files
samtools view -H ../../../data/bernardo/raw/uky_aged_stringent/minimap2/PAM54401_1271_nanopore_mapped_filtered_sorted.bam > header_1271.txt
samtools view -H ../../../data/bernardo/raw/uky_aged_stringent/minimap2/PAM54902_1291_nanopore_mapped_filtered_sorted.bam > header_1291.txt
samtools view -H ../../../data/bernardo/raw/uky_aged_stringent/minimap2/PAM54788_1304_nanopore_mapped_filtered_sorted.bam > header_1304.txt
samtools view -H ../../../data/bernardo/raw/uky_aged_stringent/minimap2/PAM54335_356_nanopore_mapped_filtered_sorted.bam > header_356.txt


## Put header back into sam files
cat header_1271.txt new_mito_reads_uky_1271.sam > new_mito_reads_uky_1271_with_header.sam
cat header_1291.txt new_mito_reads_uky_1291.sam > new_mito_reads_uky_1291_with_header.sam
cat header_1304.txt new_mito_reads_uky_1304.sam > new_mito_reads_uky_1304_with_header.sam
cat header_356.txt new_mito_reads_uky_356.sam > new_mito_reads_uky_356_with_header.sam



## Move resulting bam files to proper place
mv *header.sam ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/

## Delete intermediate files
yes | rm -r *read_ids.txt *.sam header*
