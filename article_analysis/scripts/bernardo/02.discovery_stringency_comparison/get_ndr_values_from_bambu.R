#!/usr/bin/Rscript


## Import library
library("bambu")
library("tidyverse")

## Take arguments
args <- commandArgs(trailingOnly = TRUE)
rds_name <- args[1]
output_name = args[2]

## Read RDS file
se = readRDS(rds_name)

## Only keep new transcripts
se_filtered =  se[(!is.na(mcols(se)$NDR))]

## Only keep relevant data
se_filtered_output = mcols(se_filtered)


## Drop eqClass column
drop <- c("eqClass", "eqClassById")
df_output = se_filtered_output[,!(names(se_filtered_output) %in% drop)]


# Write output
write.table(df_output, file = output_name, row.names=FALSE, sep="\t", quote = FALSE)
