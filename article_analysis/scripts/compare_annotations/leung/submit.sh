#!/bin/bash
## Submit comparison between glinos paper and our new transcripts

singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/other/novel_only_annotations/leung_annotation_clean.gtf \
    ../../../data/processed/other/novel_only_annotations/high_confidence_novel_bambu.gtf \
    -o our_high_confidence_vs_leung

