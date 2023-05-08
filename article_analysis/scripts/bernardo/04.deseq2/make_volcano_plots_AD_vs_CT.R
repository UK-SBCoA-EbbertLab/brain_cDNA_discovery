#!/usr/bin/env Rscript

## Load libraries
library(ggplot2)
library(EnhancedVolcano)

## Load datasets for differential expression results in genes with and without covariates
gene <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/genes_annotated_results_AD_vs_CT.tsv', sep="\t") 
gene_two <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/genes_AD_vs_CT_results_with_covariate_two_cell_types.tsv')
gene_four <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/genes_AD_vs_CT_results_with_covariate_four_cell_types.tsv')


## Load datasets for differential expression results in transcripts with and without covariates
transcript <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/multiple_transcripts_results_AD_vs_CT.tsv')
transcript_two <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/multiple_transcripts_results_AD_vs_CT_with_covariate_two_cell_types.tsv')
transcript_four <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/multiple_transcripts_results_AD_vs_CT_with_covariate_four_cell_types.tsv')


## Load datasets for differential expression results in transcripts with and without covariates
transcript_med <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/multiple_transcripts_results_med_relevant_AD_vs_CT.tsv')
transcript_med_two <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/multiple_transcripts_results_med_relevant_AD_vs_CT_with_covariate_two_cell_types.tsv')
transcript_med_four <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/multiple_transcripts_results_med_relevant_AD_vs_CT_with_covariate_four_cell_types.tsv')

head(gene)

head(gene$gene_name)

# PDF device
pdf("./test.pdf")

## Make fancy volcano plot and save
volcano_gene = EnhancedVolcano(gene,
                                x = 'log2FoldChange',
                                y = 'padj',
                                lab = gene$gene_name,
                                pCutoff = 0.05,
                                ylab = bquote(~-Log[10] ~ italic(Q)),
                                FCcutoff = 1,
                                pointSize = 3.0,
                                colAlpha = 0.6,
                                drawConnectors = TRUE,
                                labSize = 1.0,
                                labCol = 'black',
                                labFace = 'bold',
                                widthConnectors = 0.75,
                                colConnectors = 'black'
                                )  
                                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                plot.title = "Gene Level No Covariate AD vs CT \n",
                                plot.subtitle = element_blank(),
                                plot.caption = element_blank())


print(volcano_gene)
dev.off()
