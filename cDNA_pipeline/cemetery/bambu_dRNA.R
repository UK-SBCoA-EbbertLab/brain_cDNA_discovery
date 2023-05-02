#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

bam <- unlist(strsplit(args[1], ","))
fa_file <- args[2]
gtf_file <- args[3]

bambuAnnotations <- prepareAnnotations(gtf_file)

se_novel <- bambu(reads = bam, annotations = bambuAnnotations, genome = fa_file,
                  ncore=8, opt.discovery = list(min.sampleNumber = 2), discovery=TRUE, quant=TRUE)

writeBambuOutput(se_novel, path = "./")
