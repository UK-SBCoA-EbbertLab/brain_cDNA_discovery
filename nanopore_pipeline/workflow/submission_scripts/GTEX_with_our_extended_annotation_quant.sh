#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "../../../../../../../scratch/bag222/data/ont_data/GTEx_data_maddy/sequence_data/PCR_cDNA/*.fastq" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../../article_analysis/data/raw/nextflow_pipeline_output/bambu_discovery/extended_annotations.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS109" \
    --out_dir "./GTEX_with_our_extended_annotation_quant/" \
    --bambu_track_reads "True" \
    --is_discovery "False" \
    --is_chm13 "False" -resume
