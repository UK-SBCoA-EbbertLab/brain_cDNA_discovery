#!/bin/bash
## Submit comparison between glinos paper and our new transcripts

singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/other/novel_only_annotations/glinos_annotation_clean.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o ours_all_vs_glinos

