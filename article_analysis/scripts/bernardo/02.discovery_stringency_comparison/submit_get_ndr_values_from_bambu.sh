#!/bin/bash


singularity exec /scratch/bag222/cDNA_pipeline/singularity_container/cdna_nanopore_pipe.sif ./get_ndr_values_from_bambu.R ../../../data/bernardo/raw/merged_second_half_loose/bambu_discovery/final_discovery.RDS merged_second_half_loose_ndr.tsv
singularity exec /scratch/bag222/cDNA_pipeline/singularity_container/cdna_nanopore_pipe.sif ./get_ndr_values_from_bambu.R ../../../data/bernardo/raw/merged_first_half_loose/bambu_discovery/final_discovery.RDS merged_first_half_loose_ndr.tsv


singularity exec /scratch/bag222/cDNA_pipeline/singularity_container/cdna_nanopore_pipe.sif ./get_ndr_values_from_bambu.R ../../../data/bernardo/raw/cshl_aged_loose/bambu_discovery/final_discovery.RDS cshl_aged_loose_ndr.tsv
singularity exec /scratch/bag222/cDNA_pipeline/singularity_container/cdna_nanopore_pipe.sif ./get_ndr_values_from_bambu.R ../../../data/bernardo/raw/uky_aged_loose/bambu_discovery/final_discovery.RDS uky_aged_loose_ndr.tsv
