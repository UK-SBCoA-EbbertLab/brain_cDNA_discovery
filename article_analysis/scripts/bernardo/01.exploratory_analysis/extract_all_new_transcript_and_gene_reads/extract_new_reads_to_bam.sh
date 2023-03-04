#!/bin/bash

## Only keep the first column that contains read ids and skip the header
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_1092_PAM41667_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1092_PAM41667_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_1131_PAM44580_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1131_PAM44580_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_1163_PAM44604_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1163_PAM44604_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_1186_PAM43869_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1186_PAM43869_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_1218_PAM43779_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1218_PAM43779_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_1271_PAM44815_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1271_PAM44815_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_1291_PAG71816_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1291_PAG71816_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_1304_PAM44487_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1304_PAM44487_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_5292_PAG75292_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_5292_PAG75292_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_5295_PAG77944_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_5295_PAG77944_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_5356_PAM42933_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_5356_PAM42933_mito_read_ids.txt"}'
tail -n +2 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/sample_579_PAG75663_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "sample_579_PAG75663_mito_read_ids.txt"}'





## Extract the read ids mapping to new mitochondrial transcripts from the BAM files
samtools view ../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/QC/ | fgrep -w -f uky_1271_mito_read_ids.txt > new_mito_reads_uky_1271.sam
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
