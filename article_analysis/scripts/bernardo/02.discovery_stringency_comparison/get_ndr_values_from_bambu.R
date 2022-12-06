#!/usr/bin/Rscript


## Import library
library("bambu")

## Take arguments
args <- commandArgs(trailingOnly = TRUE)
rds_name <- args[1]
output_name = args[2]

## Read RDS file
se = readRDS(rds_name)

## Only keep new transcripts
se_filtered =  se[(!is.na(mcols(se)$txNDR))]

## Only keep relevant data
se_filtered_output = mcols(se_filtered)

# Write output
write.table(se_filtered_output, file = output_name, row.names=FALSE, sep="\t", quote = FALSE)
