#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

bam <- unlist(strsplit(args[1], ","))
fa_file <- args[2]
gtf_file <- args[3]
id <- args[4]

bambuAnnotations <- prepareAnnotations(gtf_file)

se_novel <- bambu(reads = bam, annotations = bambuAnnotations, genome = fa_file,
                  ncore=8, lowMemory = TRUE)

writeBambuOutput(se_novel, path = "./merged_final/")
