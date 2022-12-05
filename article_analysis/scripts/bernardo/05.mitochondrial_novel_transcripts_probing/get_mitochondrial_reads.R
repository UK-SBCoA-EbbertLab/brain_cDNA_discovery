#!/usr/bin/Rscript

## Import Libraries
library(bambu)
library(tidyverse)

## Read bambu expression object
se = readRDS("../../data/raw/merged_aged_stringent/bambu_discovery/final_discovery.RDS")

## Get sample names from read ids object
names_samples = names(metadata(se)$readToTranscriptMaps)

## Loop through the 4 samples
for (i in 1:4) {

        ## Extract read ids from the novel mitochondrial transcripts of interest
        read_tracking = select(metadata(se)$readToTranscriptMaps[[i]], -compatible_matches)
	read_tracking_unnested = unnest(read_tracking, equal_matches)
	read_tracking_mito = filter(read_tracking_unnested, equal_matches == 766 | equal_matches == 701 | equal_matches == 702 | equal_matches == 767)

       
        ## Write object with read ids mapping to the new mitochondrial transcripts for each sample into a csv file
        output_name = paste("./", names_samples[i], "_novel_mitochondrial_read_ids.txt", sep = "")
        print(output_name)
        write_tsv(read_tracking_mito, output_name)

}

