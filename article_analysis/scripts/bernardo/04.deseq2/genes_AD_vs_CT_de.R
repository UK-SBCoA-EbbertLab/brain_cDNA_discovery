#!/usr/bin/Rscript

## Import Libraries
library("DESeq2")



## Import counts matrix
cts <- as.matrix(read.csv("../../../data/bernardo/processed/04.deseq2/gene_counts_unfiltered.tsv", sep="\t", row.names="gene_id"))


## Import metadata
coldata <- read.csv("../../../data/bernardo/processed/04.deseq2/experimental_design.tsv", sep="\t", row.names=1)

coldata$condition <- factor(coldata$condition)
coldata$sex <- factor(coldata$sex)


## Make sure order of columns in metadata and counts matrix is the same... TRUE = good, FALSE = bad
all(rownames(coldata) == colnames(cts))



## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

## DE analysis
dds <- DESeq(dds)

## Filter genes with low counts
keep <- rowSums(counts(dds) >= 25) >= 4
dds <- dds[keep,]

## Get results
res <- results(dds)

## Summarize results
summary(res)


## Write results to tsv
write.table(res, file = "../../../data/bernardo/processed/04.deseq2/genes_AD_vs_CT_results.tsv", sep="\t")











