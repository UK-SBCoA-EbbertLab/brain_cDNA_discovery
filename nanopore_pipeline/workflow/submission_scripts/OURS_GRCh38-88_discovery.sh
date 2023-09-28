#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "../../../../../../../scratch/bag222/data/ont_data/09-02-2022_uky_6ad_6ct/*.fastq" \
    --ont_reads_txt "../../../../../../scratch/bag222/data/ont_data/09-02-2022_uky_6ad_6ct/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.88_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./OURS_GRCh38-88_discovery/" \
    --bambu_track_reads "True" \
    --is_discovery "True" \
    --is_chm13 "False"
