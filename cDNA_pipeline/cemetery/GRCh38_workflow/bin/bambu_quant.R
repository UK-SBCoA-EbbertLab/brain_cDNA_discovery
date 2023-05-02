#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

bam <- unlist(strsplit(args[1], ","))
fa_file <- args[2]
gtf_file <- args[3]
downsample <- args[4]

bambuAnnotations <- prepareAnnotations(gtf_file)

se_quant <- bambu(reads = bam, annotations = bambuAnnotations, genome = fa_file,
                  ncore=8, lowMemory = TRUE, discovery=FALSE)

output_path = paste("./bambu_quant", downsample, sep="_")

writeBambuOutput(se_quant, path = output_path)
