#!/bin/bash


singularity exec ../../../../../singularity_containers/bernardo_article_analysis.sif ./make_bambu_new_transcript_events_annotation.R \
    ../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/bambu_discovery/final_discovery.RDS \
    ../../../../data/bernardo/processed/03.gene_and_transcripts_descriptive_stats/novel_events.tsv
