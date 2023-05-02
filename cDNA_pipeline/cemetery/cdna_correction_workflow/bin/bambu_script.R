#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

bam <- unlist(strsplit(args[1], ","))
fa_file <- args[2]
gtf_file <- args[3]
id <- args[4]

print(bam)

print(bam[1])

print(bam[2])

print(typeof(bam))

print(length(bam))

bambuAnnotations <- prepareAnnotations(gtf_file)

#reads_merged <- c(bam)

se_novel <- bambu(reads = bam, annotations = bambuAnnotations, genome = fa_file,
                  ncore=8, lowMemory = TRUE)

writeBambuOutput(se_novel, path = "./merged_final/")
