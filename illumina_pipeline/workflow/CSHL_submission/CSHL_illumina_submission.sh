#!/usr/bin/env bash

## Run nextflow pipeline
nextflow ../main.nf --illumina_data "/scratch/bag222/data/technical_paper_data/cshl_illumina_data/*_R{1,2}*.fastq.gz" \
    --ref "../../../cDNA_pipeline/references/Homo_sapiens.GRCh38_ERCC.fa" \
    --transcriptome "../../../article_analysis/data/raw/nextflow_pipeline_output/transcriptome/transcriptome.fa" \
    --annotation "../../../article_analysis/data/raw/nextflow_pipeline_output/bambu_discovery/extended_annotations.gtf" \
    --housekeeping "../../../cDNA_pipeline/references/hg38.HouseKeepingGenes.bed" \
    --out_dir "./CSHL_uky_aged_brain_with_our_extended_annotation/" \
    --overhang "149" -resume

