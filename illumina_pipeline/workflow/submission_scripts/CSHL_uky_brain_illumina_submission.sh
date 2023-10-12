#!/usr/bin/env bash

## Run nextflow pipeline
nextflow ../main.nf --illumina_data "/../../../../../scratch/bag222/data/illumina_data/cshl_illumina_aged_data/*_R{1,2}*.fastq.gz" \
    --ref "../../../nanopore_pipeline/references/Homo_sapiens.GRCh38_ERCC.fa" \
    --transcriptome "../../../article_analysis/data/raw/nextflow_pipeline_output/transcriptome/transcriptome.fa" \
    --annotation "../../../article_analysis/data/raw/nextflow_pipeline_output/bambu_discovery/extended_annotations.gtf" \
    --housekeeping "../../../cDNA_pipeline/references/hg38.HouseKeepingGenes.bed" \
    --out_dir "./CSHL_illumina_uky_aged_brain_with_our_extended_annotation/" \
    --overhang "149" -resume

