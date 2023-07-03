#!/usr/bin/Rscript

## Import Libraries
library(bambu)
library(tidyverse)

## Read bambu expression object
se = readRDS("../../../../data/bernardo/raw/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/bambu_discovery/final_discovery.RDS")

## Get sample names from read ids object
names_samples = names(metadata(se)$readToTranscriptMaps)

## Loop through the 12 samples
for (i in 1:12) {

        ## Extract read ids from the novel mitochondrial transcripts of interest
        read_tracking = select(metadata(se)$readToTranscriptMaps[[i]], -compatibleMatches)
	read_tracking_unnested = unnest(read_tracking, equalMatches)
	read_tracking_mito = filter(read_tracking_unnested, equalMatches < 3430)

        ## Write object with read ids mapping to the new mitochondrial transcripts for each sample into a csv file
        output_name = paste("./", names_samples[i], "_novel_read_ids.txt", sep = "")
        print(output_name)
        write_tsv(read_tracking_mito, output_name)

}

