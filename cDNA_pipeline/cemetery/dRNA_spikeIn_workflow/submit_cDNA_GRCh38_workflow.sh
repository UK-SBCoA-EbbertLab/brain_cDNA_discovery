#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_test_direct_rna_data/*ERCC*.fastq" \
    --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_test_direct_rna_data/*ERCC*.txt" \
    --ref "../../references/ERCC92.fa" \
    --annotation "../../references/ERCC92.gtf" \
    --out_dir "./spikeIn/" \
    --is_chm13 "False" \
    --is_dRNA "True" \
    --is_spikeIn "True" -resume
