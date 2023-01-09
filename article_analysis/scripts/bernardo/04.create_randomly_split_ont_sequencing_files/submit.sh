#!/bin/bash


sbatch ./split_flowcell.ogs /scratch/bag222/data/ont_data/uky_aged_data/PAM54335_356_nanopore.txt \
    /scratch/bag222/data/ont_data/uky_aged_data/PAM54335_356_nanopore.fastq


sbatch ./split_flowcell.ogs /scratch/bag222/data/ont_data/uky_aged_data/PAM54401_1271_nanopore.txt \
    /scratch/bag222/data/ont_data/uky_aged_data/PAM54401_1271_nanopore.fastq


sbatch ./split_flowcell.ogs /scratch/bag222/data/ont_data/uky_aged_data/PAM54788_1304_nanopore.txt \
    /scratch/bag222/data/ont_data/uky_aged_data/PAM54788_1304_nanopore.fastq


sbatch ./split_flowcell.ogs /scratch/bag222/data/ont_data/uky_aged_data/PAM54902_1291_nanopore.txt \
    /scratch/bag222/data/ont_data/uky_aged_data/PAM54902_1291_nanopore.fastq

