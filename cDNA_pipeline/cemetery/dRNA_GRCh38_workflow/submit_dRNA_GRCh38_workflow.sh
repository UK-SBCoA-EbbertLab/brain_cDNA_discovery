#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_test_direct_rna_data/*.fastq" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.104.gtf" \
    --out_dir "./GRCh38_dRNA/" \
    --is_chm13 "False" \
    --is_dRNA "True" -resume
