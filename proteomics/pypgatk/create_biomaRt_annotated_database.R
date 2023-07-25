library(biomaRt)
library(tidyverse)

f <- file("stdin")

transcript_ids<-read.table(f, header = FALSE, sep = "\t")

#useful information about biomaRt: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
# and more: https://bioc.ism.ac.jp/packages/3.4/bioc/manuals/biomaRt/man/biomaRt.pdf

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 107)

protein_seqs<-getSequence(id=transcript_ids[,1], type = "ensembl_transcript_id", seqType = "peptide", mart = ensembl)

exportFASTA(protein_seqs, "")

