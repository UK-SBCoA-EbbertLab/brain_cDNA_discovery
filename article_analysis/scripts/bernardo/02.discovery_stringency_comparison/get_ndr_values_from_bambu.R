#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

rds_name <- args[1]
output_name = args[2]

se = readRDS(rds_name)

se_filtered =  se[(!is.na(mcols(se)$txNDR))]

se_filtered_output = mcols(se_filtered)

write.table(se_filtered_output, file = output_name, row.names=FALSE, sep="\t", quote = FALSE)
