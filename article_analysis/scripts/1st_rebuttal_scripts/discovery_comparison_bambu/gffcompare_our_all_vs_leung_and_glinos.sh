#!/bin/bash

## OURS ALL vs GTEx ENSEMBL 107 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/other/novel_only_annotations/leung_annotation_clean.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o OURS_ALL_vs_leung

## OURS ALL vs GTEx ENSEMBL 88 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/other/novel_only_annotations/glinos_annotation_clean.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o OURS_ALL_vs_glinos
