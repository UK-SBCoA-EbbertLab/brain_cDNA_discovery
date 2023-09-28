#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "../../../../../../../scratch/bag222/data/ont_data/09-02-2022_uky_6ad_6ct/*.fastq" \
    --ont_reads_txt "../../../../../../scratch/bag222/data/ont_data/09-02-2022_uky_6ad_6ct/*.txt" \
    --ref "../../references/chm13v2.0_ERCC.fa" \
    --annotation "../../references/CHM13.v2.0.gff3" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./OURS_CHM13_quant/" \
    --bambu_track_reads "True" \
    --is_discovery "False" \
    --ercc "../../references/ERCC92.gtf" \
    --is_chm13 "True"
