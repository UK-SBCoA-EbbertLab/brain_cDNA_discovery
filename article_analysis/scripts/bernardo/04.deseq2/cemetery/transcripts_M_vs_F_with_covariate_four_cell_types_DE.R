#!/usr/bin/Rscript

## Import Libraries
library("DESeq2")



## Import counts matrix
cts <- as.matrix(read.csv("../../../data/bernardo/processed/04.deseq2/transcript_counts_multiple_median_cpm_1.tsv", sep="\t", row.names="transcript_id"))
cts_med <- as.matrix(read.csv("../../../data/bernardo/processed/04.deseq2/transcript_counts_multiple_med_relevant_median_cpm_1.tsv", sep="\t", row.names="transcript_id"))


## Import metadata
coldata <- read.csv("../../../data/bernardo/processed/04.deseq2/experimental_design_with_four_cell_populations.tsv", sep="\t", row.names=1)

head(coldata, 12)

## Convert proportion_neuronal from numerical to factor
coldata$proportion_neuronal <- cut(coldata$proportion_neuronal, 3, labels=c('low_neu', 'average_neu', 'high_neu'))
coldata$proportion_astro <- cut(coldata$proportion_astro, 3, labels=c('low_astro', 'average_astro', 'high_astro'))
coldata$proportion_micro <- cut(coldata$proportion_micro, 3, labels=c('low_micro', 'average_micro', 'high_micro'))
coldata$proportion_oligo <- cut(coldata$proportion_oligo, 3, labels=c('low_oligo', 'average_oligo', 'high_oligo'))


coldata$condition <- factor(coldata$condition)
coldata$sex <- factor(coldata$sex)
coldata$proportion_neuronal <- factor(coldata$proportion_neuronal)
coldata$proportion_astro <- factor(coldata$proportion_astro)
coldata$proportion_micro <- factor(coldata$proportion_micro)
coldata$proportion_oligo <- factor(coldata$proportion_oligo)



## Make sure order of columns in metadata and counts matrix is the same... TRUE = good, FALSE = bad
all(rownames(coldata) == colnames(cts))
all(rownames(coldata) == colnames(cts_med))



## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ proportion_astro +  proportion_neuronal + proportion_micro + sex)

## Create DESeq2 object
dds_med <- DESeqDataSetFromMatrix(countData = cts_med,
                              colData = coldata,
                              design = ~ proportion_astro +  proportion_neuronal + proportion_micro + sex)


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
write.table(res, file = "../../../data/bernardo/processed/04.deseq2/multiple_transcripts_results_M_vs_F_with_covariate_four_cell_types.tsv", sep="\t")
write.table(res_med, file =  "../../../data/bernardo/processed/04.deseq2/multiple_transcripts_results_med_relevant_M_vs_F_with_covariate_four_cell_types.tsv", sep="\t")











