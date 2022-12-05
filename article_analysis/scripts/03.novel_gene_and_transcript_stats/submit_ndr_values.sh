#!/bin/bash


singularity exec /scratch/bag222/cDNA_pipeline/singularity_container/cdna_nanopore_pipe.sif ./ndr_values.R ../../data/raw/merged_aged_stringent/bambu_discovery/final_discovery.RDS merged_aged_stringent_ndr.tsv
