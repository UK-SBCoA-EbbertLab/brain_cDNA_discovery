#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "../../../../../../../scratch/bag222/data/ont_data/GTEx_data_maddy/sequence_data/PCR_cDNA/*.fastq" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.88_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS109" \
    --out_dir "./GTEX_with_GRCh38-88_discovery/" \
    --bambu_track_reads "True" \
    --is_discovery "True" \
    --is_chm13 "False" -resume
