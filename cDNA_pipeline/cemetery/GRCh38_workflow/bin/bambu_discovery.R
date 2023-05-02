#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

bam <- unlist(strsplit(args[1], ","))
fa_file <- args[2]
gtf_file <- args[3]

bambuAnnotations <- prepareAnnotations(gtf_file)

help(bambu, package="bambu")

se_novel <- bambu(reads = bam, annotations = bambuAnnotations, genome = fa_file,
                  ncore=8, opt.discovery = list(min.sampleNumber = 2))

writeBambuOutput(se_novel, path = "./bambu_discovery/")
