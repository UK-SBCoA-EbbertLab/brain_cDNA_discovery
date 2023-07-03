#!/usr/bin/env Rscript

## Load bambu library
library(bambu)
library(tidyverse)
library(dplyr)

## Get arguments from command line
args <- commandArgs(trailingOnly = TRUE)

input_name <- args[1]
output_name <- args[2]


## Load bambu R object
se = readRDS(input_name)

## Create dataframe containing transcript descriptions
se_novel_only = se[rowData(se)$txClassDescription != "annotation", ]
se_novel_only_2 = rowData(se_novel_only)

## Select relevant columns for output
se_novel_only_final = se_novel_only_2[c("GENEID", "TXNAME", "txClassDescription")]

## Convert to dataframe
df_output <- as_tibble(se_novel_only_final)

## Write tsv with output
write_tsv(df_output, output_name)
