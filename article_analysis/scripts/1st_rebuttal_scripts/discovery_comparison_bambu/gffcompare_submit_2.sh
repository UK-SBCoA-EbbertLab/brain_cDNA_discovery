#!/bin/bash

## GTEX ENSEMBL 88 vs GTEx ENSEMBL 107 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_88_NEW_ALL.gtf \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_107_NEW_ALL.gtf \
    -o GTEx_ENSEMBL_88_ALL_vs_GTEx_ENSEMBL_107_ALL

## OURS ENSEMBL 88 vs OURS ENSEMBL 107 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_88_NEW_ALL.gtf \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_107_NEW_ALL.gtf \
    -o OURS_ENSEMBL_88_ALL_vs_OURS_ENSEMBL_107_ALL


## GTEX ENSEMBL 107 vs GTEx ENSEMBL 88 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_107_NEW_ALL.gtf \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_88_NEW_ALL.gtf \
    -o GTEx_ENSEMBL_107_ALL_vs_GTEx_ENSEMBL_88_ALL

## OURS ENSEMBL 107 vs OURS ENSEMBL 88 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_107_NEW_ALL.gtf \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_88_NEW_ALL.gtf \
    -o OURS_ENSEMBL_107_ALL_vs_OURS_ENSEMBL_88_ALL
