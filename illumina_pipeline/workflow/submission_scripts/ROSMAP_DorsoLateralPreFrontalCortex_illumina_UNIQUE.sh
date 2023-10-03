#!/usr/bin/env bash

## Run nextflow pipeline
nextflow ../main.nf --bam "/../../../../../scratch/bag222/data/illumina_data/ROSMAP_illumina_DorsoLateralPreFrontalCortex_with_our_extended_annotation_BAMS/RISK_255_S78_Aligned.toTranscriptome.out.bam" \
    --transcriptome "../../../article_analysis/data/raw/nextflow_pipeline_output/transcriptome/transcriptome.fa" \
    --out_dir "./ROSMAP_illumina_DorsoLateralPreFrontalCortex_UNIQUE/" -resume

