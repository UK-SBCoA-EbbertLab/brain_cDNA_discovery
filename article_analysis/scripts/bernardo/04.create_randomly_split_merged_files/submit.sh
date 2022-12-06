#!/bin/bash


sbatch ./split_flowcell.ogs /scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data/cshl_1271_uky.txt \
    /scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data/cshl_1271_uky.fastq


sbatch ./split_flowcell.ogs /scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data/cshl_1291_uky.txt \
    /scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data/cshl_1291_uky.fastq


sbatch ./split_flowcell.ogs /scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data/cshl_1304_uky.txt \
    /scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data/cshl_1304_uky.fastq


sbatch ./split_flowcell.ogs /scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data/cshl_356_uky.txt \
    /scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data/cshl_356_uky.fastq

