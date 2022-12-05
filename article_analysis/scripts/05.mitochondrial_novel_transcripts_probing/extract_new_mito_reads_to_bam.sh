#!/bin/bash


tail -n +2 ../../data/processed/05.mitochondrial_novel_transcripts_probing/cshl_1271_uky_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "cshl_1271_uky_mito_read_ids.txt"}'
tail -n +2 ../../data/processed/05.mitochondrial_novel_transcripts_probing/cshl_1291_uky_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "cshl_1291_uky_mito_read_ids.txt"}'
tail -n +2 ../../data/processed/05.mitochondrial_novel_transcripts_probing/cshl_1304_uky_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "cshl_1304_uky_mito_read_ids.txt"}'
tail -n +2 ../../data/processed/05.mitochondrial_novel_transcripts_probing/cshl_356_uky_mapped_filtered_sorted_novel_mitochondrial_read_ids.txt | awk  -F'\t' '{print $1 > "cshl_356_uky_mito_read_ids.txt"}'


samtools view <PATH TO 1271 BAM> | fgrep -w -f cshl_1271_uky_mito_read_ids.txt > new_mito_reads_cshl_1271_uky.bam
samtools view <PATH TO 1291 BAM> | fgrep -w -f cshl_1291_uky_mito_read_ids.txt > new_mito_reads_cshl_1291_uky.bam
samtools view <PATH TO 1304 BAM> | fgrep -w -f cshl_1304_uky_mito_read_ids.txt > new_mito_reads_cshl_1304_uky.bam
samtools view <PATH TO 356 BAM> | fgrep -w -f cshl_356_uky_mito_read_ids.txt > new_mito_reads_cshl_356_uky.bam


mv *.bam ../../data/processed/05.mitochondrial_novel_transcripts_probing/

yes | rm -r *read_ids.txt
