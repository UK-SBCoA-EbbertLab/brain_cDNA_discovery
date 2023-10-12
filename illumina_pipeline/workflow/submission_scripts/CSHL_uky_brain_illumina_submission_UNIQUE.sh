#!/usr/bin/env bash

## Run nextflow pipeline
nextflow ../main.nf --bam "./results/CSHL_illumina_uky_aged_brain_with_our_extended_annotation/STAR/*.toTranscriptome.out.bam" \
    --transcriptome "../../../article_analysis/data/raw/nextflow_pipeline_output/transcriptome/transcriptome.fa" \
    --out_dir "./CSHL_illumina_uky_aged_brain_with_our_extended_annotation_UNIQUE/" -resume

