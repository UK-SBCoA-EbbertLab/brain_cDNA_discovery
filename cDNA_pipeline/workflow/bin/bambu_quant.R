#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

rc_files <- unlist(strsplit(args[1], ","))
fa_file <- args[2]
gtf_file <- args[3]

bambuAnnotations <- prepareAnnotations(gtf_file)

se_quant <- bambu(reads=rc_files, annotations=bambuAnnotations, genome=fa_file,
                  ncore=1, lowMemory=TRUE, discovery=FALSE, quant=TRUE)

writeBambuOutput(se_quant, path = "./bambu_quant/")

saveRDS(se_quant, file="./bambu_quant/final_quant.RDS")
