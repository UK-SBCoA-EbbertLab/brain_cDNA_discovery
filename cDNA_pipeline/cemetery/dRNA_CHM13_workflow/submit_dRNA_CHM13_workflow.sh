#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_test_direct_rna_data/*.fastq" \
    --ref "../../references/chm13v2.0.fa" \
    --annotation "../../references/CHM13.v2.0.gff3" \
    --out_dir "./CHM13_dRNA/" \
    --is_chm13 "True" \
    --is_dRNA "True" -resume
