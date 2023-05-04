#!/usr/bin/Rscript

## Import Libraries
library("DESeq2")



## Import counts matrix
cts <- as.matrix(read.csv("../../../data/bernardo/processed/04.deseq2/gene_counts_unfiltered.tsv", sep="\t", row.names="gene_id"))


## Import metadata
coldata <- read.csv("../../../data/bernardo/processed/04.deseq2/experimental_design_with_two_cell_populations.tsv", sep="\t", row.names=1)

head(coldata, 12)

## Convert proportion_neuronal from numerical to factor
coldata$proportion_neuronal <- cut(coldata$proportion_neuronal, 3, labels=c('neu_low', 'neu_avg', 'new_high'))
coldata$proportion_non_neuronal <- cut(coldata$proportion_non_neuronal, 3, labels=c('low_other', 'average_other', 'high_other'))


coldata$condition <- factor(coldata$condition)
coldata$sex <- factor(coldata$sex)
coldata$proportion_neuronal <- factor(coldata$proportion_neuronal)
coldata$proportion_non_neuronal <- factor(coldata$proportion_non_neuronal)

head(coldata, 12)

## Make sure order of columns in metadata and counts matrix is the same... TRUE = good, FALSE = bad
all(rownames(coldata) == colnames(cts))



## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ proportion_neuronal + condition)

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
write.table(res, file = "../../../data/bernardo/processed/04.deseq2/genes_AD_vs_CT_results_with_covariate_two_cell_types.tsv", sep="\t")











