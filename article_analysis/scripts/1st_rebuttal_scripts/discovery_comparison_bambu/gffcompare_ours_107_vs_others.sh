#!/bin/bash

## OURS ALL vs GTEx ENSEMBL 107 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_107_NEW_ALL.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o OURS_ALL_vs_GTEx_ENSEMBL_107_ALL

## OURS ALL vs GTEx ENSEMBL 88 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_88_NEW_ALL.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o OURS_ALL_vs_GTEx_ENSEMBL_88_ALL

## OURS ALL vs GTEx ENSEMBL 107 HF
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_107_NEW_HF.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o OURS_ALL_vs_GTEx_ENSEMBL_107_HF

## OURS ALL vs GTEx ENSEMBL 88 HF
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_88_NEW_HF.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o OURS_ALL_vs_GTEx_ENSEMBL_88_HF

## OURS ALL vs OURS ENSEMBL 107 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_107_NEW_ALL.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o OURS_ALL_vs_OURS_ENSEMBL_107_ALL

## OURS ALL vs OURS ENSEMBL 88 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_88_NEW_ALL.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o OURS_ALL_vs_OURS_ENSEMBL_88_ALL

## OURS ALL vs OURS ENSEMBL 107 HF
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_107_NEW_HF.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o OURS_ALL_vs_OURS_ENSEMBL_107_HF

## OURS ALL vs OURS ENSEMBL 88 HF
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_88_NEW_HF.gtf \
    ../../../data/processed/other/novel_only_annotations/all_novel_bambu.gtf \
    -o OURS_ALL_vs_OURS_ENSEMBL_88_HF






## OURS HF vs GTEx ENSEMBL 107 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_107_NEW_ALL.gtf \
    ../../../data/processed/other/novel_only_annotations/high_confidence_novel_bambu.gtf \
    -o OURS_HF_vs_GTEx_ENSEMBL_107_ALL

## OURS HF vs GTEx ENSEMBL 88 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_88_NEW_ALL.gtf \
    ../../../data/processed/other/novel_only_annotations/high_confidence_novel_bambu.gtf \
    -o OURS_HF_vs_GTEx_ENSEMBL_88_ALL

## OURS HF vs GTEx ENSEMBL 107 HF
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_107_NEW_HF.gtf \
    ../../../data/processed/other/novel_only_annotations/high_confidence_novel_bambu.gtf \
    -o OURS_HF_vs_GTEx_ENSEMBL_107_HF

## OURS HF vs GTEx ENSEMBL 88 HF
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/GTEx_ENSEMBL_88_NEW_HF.gtf \
    ../../../data/processed/other/novel_only_annotations/high_confidence_novel_bambu.gtf \
    -o OURS_HF_vs_GTEx_ENSEMBL_88_HF

## OURS HF vs OURS ENSEMBL 107 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_107_NEW_ALL.gtf \
    ../../../data/processed/other/novel_only_annotations/high_confidence_novel_bambu.gtf \
    -o OURS_HF_vs_OURS_ENSEMBL_107_ALL

## OURS HF vs OURS ENSEMBL 88 ALL
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_88_NEW_ALL.gtf \
    ../../../data/processed/other/novel_only_annotations/high_confidence_novel_bambu.gtf \
    -o OURS_HF_vs_OURS_ENSEMBL_88_ALL

## OURS HF vs OURS ENSEMBL 107 HF
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_107_NEW_HF.gtf \
    ../../../data/processed/other/novel_only_annotations/high_confidence_novel_bambu.gtf \
    -o OURS_HF_vs_OURS_ENSEMBL_107_HF

## OURS HF vs OURS ENSEMBL 88 HF
singularity exec ../../../../singularity_containers/nanopore.sif gffcompare -r \
    ../../../data/processed/1st_rebuttal/discovery_comparison_bambu/OURS_ENSEMBL_88_NEW_HF.gtf \
    ../../../data/processed/other/novel_only_annotations/high_confidence_novel_bambu.gtf \
    -o OURS_HF_vs_OURS_ENSEMBL_88_HF



