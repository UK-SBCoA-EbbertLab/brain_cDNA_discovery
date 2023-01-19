#!/bin/bash


singularity exec ../../../../singularity_containers/bambu.sif ./get_ndr_values_from_bambu.R ../../../data/bernardo/raw/uky_aged_firstHalf_loose/bambu_discovery/final_discovery.RDS uky_firstHalf_loose_ndr.tsv
singularity exec ../../../../singularity_containers/bambu.sif ./get_ndr_values_from_bambu.R ../../../data/bernardo/raw/uky_aged_secondHalf_loose/bambu_discovery/final_discovery.RDS uky_secondHalf_loose_ndr.tsv
