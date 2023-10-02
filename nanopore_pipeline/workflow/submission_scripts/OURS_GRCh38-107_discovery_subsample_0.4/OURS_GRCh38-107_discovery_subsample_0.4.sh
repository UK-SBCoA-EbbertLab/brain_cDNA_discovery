#!/usr/bin/env bash


nextflow ../../main.nf --ont_reads_fq "../../../../../../../../scratch/bag222/data/ont_data/make_downsampled_fastqs/data/*-0.4.fastq" \
    --ont_reads_txt "../../../../../../../../scratch/bag222/data/ont_data/make_downsampled_fastqs/data/*-0.4.txt" \
    --ref "../../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../../references/hg38.HouseKeepingGenes.bed" \
    --multiqc_config "../../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./OURS_GRCh38-107_discovery_subsample_0.4/" \
    --bambu_track_reads "True" \
    --is_discovery "True" \
    --is_chm13 "False" -resume
