#!/bin/bash


## Script takes: gtf for comparison, gtf to be compared, name of altenative column (not novel), output file name

## Output file will contain all transcripts from the gtf to be compared, including whether they are novel or annotated


## Loose run comparison, bambu NDR = 1
sbatch run_identify_novel_transcripts.ogs ../../../data/bernardo/processed/99.other/create_annotations/novel_only_annotations/uky_aged_firstHalf_loose_annotation_novel_only.gtf \
    ../../../data/bernardo/processed/99.other/create_annotations/novel_only_annotations/uky_aged_secondHalf_loose_annotation_novel_only.gtf "ANNOTATED" second_half_annotation_loose.tsv

sbatch run_identify_novel_transcripts.ogs ../../../data/bernardo/processed/99.other/create_annotations/novel_only_annotations/uky_aged_secondHalf_loose_annotation_novel_only.gtf \
    ../../../data/bernardo/processed/99.other/create_annotations/novel_only_annotations/uky_aged_firstHalf_loose_annotation_novel_only.gtf "ANNOTATED" first_half_annotation_loose.tsv
