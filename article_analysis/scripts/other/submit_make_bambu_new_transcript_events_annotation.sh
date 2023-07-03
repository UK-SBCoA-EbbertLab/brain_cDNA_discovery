#!/bin/bash


singularity exec ../../../singularity_containers/bernardo_article_analysis.sif ./make_bambu_new_transcript_events_annotation.R \
    ../../data/raw/nextflow_pipeline_output/bambu_discovery/final_discovery.RDS \
    ../../data/processed/paper_figures/novel_events.tsv
