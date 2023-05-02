#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/*PA*.fastq" \
    --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/*PA*.txt" \
    --ref "../../references/ERCC92.fa" \
    --annotation "../../references/ERCC92.gtf" \
    --out_dir "./spikeIn/" \
    --is_chm13 "False" \
    --is_dRNA "False" \
    --is_spikeIn "True" -resume
