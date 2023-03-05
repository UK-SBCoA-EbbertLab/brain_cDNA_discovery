#!/bin/bash

## Only keep the first column that contains read ids and skip the header
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_1092_PAM41667_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1092_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_1131_PAM44580_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1131_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_1163_PAM44604_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1163_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_1186_PAM43869_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1186_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_1218_PAM43779_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1218_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_1271_PAM44815_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1271_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_1291_PAG71816_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1291_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_1304_PAM44487_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_1304_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_5292_PAG75292_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_5292_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_5295_PAG77944_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_5295_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_5356_PAM42933_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_5356_new_read_ids.txt"}'
tail -n +2 ../../../../data/bernardo/processed/01.exploratory_analysis/sample_579_PAG75663_mapped_filtered_sorted_novel_read_ids.txt | awk  -F'\t' '{print $1 > "sample_579_new_read_ids.txt"}'



## Extract the read ids mapping to new mitochondrial transcripts from the BAM files
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1092_PAM41667_mapped_filtered_sorted.bam | fgrep -w -f sample_1092_new_read_ids.txt > new_reads_sample_1092.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1131_PAM44580_mapped_filtered_sorted.bam | fgrep -w -f sample_1131_new_read_ids.txt > new_reads_sample_1131.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1163_PAM44604_mapped_filtered_sorted.bam | fgrep -w -f sample_1163_new_read_ids.txt > new_reads_sample_1163.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1186_PAM43869_mapped_filtered_sorted.bam | fgrep -w -f sample_1186_new_read_ids.txt > new_reads_sample_1186.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1218_PAM43779_mapped_filtered_sorted.bam | fgrep -w -f sample_1218_new_read_ids.txt > new_reads_sample_1218.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1271_PAM44815_mapped_filtered_sorted.bam | fgrep -w -f sample_1271_new_read_ids.txt > new_reads_sample_1271.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1291_PAG71816_mapped_filtered_sorted.bam | fgrep -w -f sample_1291_new_read_ids.txt > new_reads_sample_1291.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1304_PAM44487_mapped_filtered_sorted.bam | fgrep -w -f sample_1304_new_read_ids.txt > new_reads_sample_1304.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_5292_PAG75292_mapped_filtered_sorted.bam | fgrep -w -f sample_5292_new_read_ids.txt > new_reads_sample_5292.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_5295_PAG77944_mapped_filtered_sorted.bam | fgrep -w -f sample_5295_new_read_ids.txt > new_reads_sample_5295.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_5356_PAM42933_mapped_filtered_sorted.bam | fgrep -w -f sample_5356_new_read_ids.txt > new_reads_sample_5356.sam
samtools view ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_579_PAG75663_mapped_filtered_sorted.bam | fgrep -w -f sample_579_new_read_ids.txt > new_reads_sample_579.sam

## Get header from bam files
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1092_PAM41667_mapped_filtered_sorted.bam > header_1092.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1131_PAM44580_mapped_filtered_sorted.bam > header_1131.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1163_PAM44604_mapped_filtered_sorted.bam > header_1163.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1186_PAM43869_mapped_filtered_sorted.bam > header_1186.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1218_PAM43779_mapped_filtered_sorted.bam > header_1218.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1271_PAM44815_mapped_filtered_sorted.bam > header_1271.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1291_PAG71816_mapped_filtered_sorted.bam > header_1291.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_1304_PAM44487_mapped_filtered_sorted.bam > header_1304.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_5292_PAG75292_mapped_filtered_sorted.bam > header_5292.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_5295_PAG77944_mapped_filtered_sorted.bam > header_5295.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_5356_PAM42933_mapped_filtered_sorted.bam > header_5356.txt
samtools view -H ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/mapping_cDNA/sample_579_PAG75663_mapped_filtered_sorted.bam > header_579.txt



## Put header back into sam files
cat header_1092.txt new_reads_sample_1092.sam > new_reads_sample_1092_with_header.sam
cat header_1131.txt new_reads_sample_1131.sam > new_reads_sample_1131_with_header.sam
cat header_1163.txt new_reads_sample_1163.sam > new_reads_sample_1163_with_header.sam
cat header_1186.txt new_reads_sample_1186.sam > new_reads_sample_1186_with_header.sam
cat header_1218.txt new_reads_sample_1218.sam > new_reads_sample_1218_with_header.sam
cat header_1271.txt new_reads_sample_1271.sam > new_reads_sample_1271_with_header.sam
cat header_1291.txt new_reads_sample_1291.sam > new_reads_sample_1291_with_header.sam
cat header_1304.txt new_reads_sample_1304.sam > new_reads_sample_1304_with_header.sam
cat header_5292.txt new_reads_sample_5292.sam > new_reads_sample_5292_with_header.sam
cat header_5295.txt new_reads_sample_5295.sam > new_reads_sample_5295_with_header.sam
cat header_5356.txt new_reads_sample_5356.sam > new_reads_sample_5356_with_header.sam
cat header_579.txt new_reads_sample_579.sam > new_reads_sample_579_with_header.sam


## Move resulting bam files to proper place
mv *header.sam ../../../../data/bernardo/processed/01.exploratory_analysis/


## Delete intermediate files
yes | rm -r *read_ids.txt *.sam header*
