#!/bin/bash
## Submit comparison between glinos paper and our new transcripts

singularity exec ../../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../../data/bernardo/processed/99.other/create_annotations/novel_only_annotations/glinos_annotation_clean.gtf \
    ../../../../data/bernardo/processed/99.other/create_annotations/novel_only_annotations/high_confidence_novel_bambu.gtf \
    -o ours_high_confidence_vs_glinos

