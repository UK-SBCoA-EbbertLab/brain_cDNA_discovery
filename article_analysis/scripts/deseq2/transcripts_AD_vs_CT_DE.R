#!/usr/bin/Rscript

## Import Libraries
library("DESeq2")



## Import counts matrix
cts <- as.matrix(read.csv("../../data/processed/deseq2/transcript_counts_multiple_median_cpm_1.tsv", sep="\t", row.names="transcript_id"))
cts_med <- as.matrix(read.csv("../../data/processed/deseq2/transcript_counts_multiple_med_relevant_median_cpm_1.tsv", sep="\t", row.names="transcript_id"))


## Import metadata
coldata <- read.csv("../../data/processed/deseq2/experimental_design.tsv", sep="\t", row.names=1)

coldata$condition <- factor(coldata$condition)
coldata$sex <- factor(coldata$sex)


## Make sure order of columns in metadata and counts matrix is the same... TRUE = good, FALSE = bad
all(rownames(coldata) == colnames(cts))
all(rownames(coldata) == colnames(cts_med))



## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

## Create DESeq2 object
dds_med <- DESeqDataSetFromMatrix(countData = cts_med,
                              colData = coldata,
                              design = ~ condition)


## DE analysis
dds <- DESeq(dds)
dds_med <- DESeq(dds_med)

## Get results
res <- results(dds)
res_med <- results(dds_med)

## Summarize results
summary(res)
summary(res_med)


## Write results to tsv
write.table(res, file = "../../data/processed/deseq2/multiple_transcripts_results_AD_vs_CT.tsv", sep="\t")
write.table(res_med, file =  "../../data/processed/deseq2/multiple_transcripts_results_med_relevant_AD_vs_CT.tsv", sep="\t")











