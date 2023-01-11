#!/bin/bash


singularity exec ../../../../../singularity_containers/bernardo_article_analysis.sif ./make_bambu_new_transcript_events_annotation.R \
    ../../../../data/bernardo/raw/uky_aged_stringent/bambu_discovery/final_discovery.RDS ../../../../data/bernardo/processed/03.gene_and_transcripts_descriptive_stats/uky_aged_stringent_novel_events.tsv
