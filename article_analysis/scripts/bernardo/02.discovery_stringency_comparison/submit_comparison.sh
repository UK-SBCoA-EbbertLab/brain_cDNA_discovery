#!/bin/bash


## Script takes: gtf for comparison, gtf to be compared, name of altenative column (not novel), output file name

## Output file will contain all transcripts from the gtf to be compared, including whether they are novel or annotated


## Loose run comparison, bambu NDR = 1
sbatch run_identify_novel_transcripts.ogs ../../data/processed/novel_only_annotations/merged_first_half_loose_annotation_novel_only.gtf \
    ../../data/processed/novel_only_annotations/merged_second_half_loose_annotation_novel_only.gtf "ANNOTATED" second_half_annotation_loose.tsv

sbatch run_identify_novel_transcripts.ogs ../../data/processed/novel_only_annotations/merged_second_half_loose_annotation_novel_only.gtf \
    ../../data/processed/novel_only_annotations/merged_first_half_loose_annotation_novel_only.gtf "ANNOTATED" first_half_annotation_loose.tsv


## Moderate run comparison, bambu NDR = 0.5
sbatch run_identify_novel_transcripts.ogs ../../data/processed/novel_only_annotations/merged_first_half_moderate_annotation_novel_only.gtf \
    ../../data/processed/novel_only_annotations/merged_second_half_moderate_annotation_novel_only.gtf "ANNOTATED" second_half_annotation_moderate.tsv

sbatch run_identify_novel_transcripts.ogs ../../data/processed/novel_only_annotations/merged_second_half_moderate_annotation_novel_only.gtf \
    ../../data/processed/novel_only_annotations/merged_first_half_moderate_annotation_novel_only.gtf "ANNOTATED" first_half_annotation_moderate.tsv


## Stringent run comparison, bambu NDR = 0.1
sbatch run_identify_novel_transcripts.ogs ../../data/processed/novel_only_annotations/merged_first_half_stringent_annotation_novel_only.gtf \
    ../../data/processed/novel_only_annotations/merged_second_half_stringent_annotation_novel_only.gtf "ANNOTATED" second_half_annotation_stringent.tsv

sbatch run_identify_novel_transcripts.ogs ../../data/processed/novel_only_annotations/merged_second_half_stringent_annotation_novel_only.gtf \
    ../../data/processed/novel_only_annotations/merged_first_half_stringent_annotation_novel_only.gtf "ANNOTATED" first_half_annotation_stringent.tsv


## Super stringent run comparison, bambu NDR = 0.01
sbatch run_identify_novel_transcripts.ogs ../../data/processed/novel_only_annotations/merged_first_half_super_stringent_annotation_novel_only.gtf \
    ../../data/processed/novel_only_annotations/merged_second_half_super_stringent_annotation_novel_only.gtf "ANNOTATED" second_half_annotation_super_stringent.tsv

sbatch run_identify_novel_transcripts.ogs ../../data/processed/novel_only_annotations/merged_second_half_super_stringent_annotation_novel_only.gtf \
    ../../data/processed/novel_only_annotations/merged_first_half_super_stringent_annotation_novel_only.gtf "ANNOTATED" first_half_annotation_super_stringent.tsv

