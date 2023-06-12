#!/usr/bin/env Rscript

## Load libraries
library(ggplot2)
library(EnhancedVolcano)

## Load datasets for differential expression results in genes with and without covariates
gene <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/genes_annotated_results_M_vs_F.tsv', sep="\t") 


## Load datasets for differential expression results in transcripts with and without covariates
transcript <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/annotated_multiple_transcripts_results_M_vs_F.tsv', sep="\t")


# PDF device
pdf("../../../figures/bernardo/04.deseq2/volcano_plot_SEX.pdf")

## Make fancy volcano plot and save
volcano_gene = EnhancedVolcano(gene,
                                x = 'log2FoldChange',
                                y = 'padj',
                                lab = gene$gene_name,
                                selectLab = c('TGOLN2', 'DYRK1B', 'CLU', 'NDUFA10'),
                                pCutoff = 0.05,
                                ylab = bquote(~-Log[10] ~ italic(Q)),
                                FCcutoff = 1,
                                pointSize = 3.0,
                                colAlpha = 0.6,
                                drawConnectors = TRUE,
                                labSize = 4.5,
                                labCol = 'black',
                                labFace = 'bold',
                                widthConnectors = 0.5,
                                boxedLabels = TRUE,
                                lengthConnectors = unit(0.01, 'npc'),
                                colConnectors = 'black'
                                ) +     
                                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                plot.title = element_blank(),
                                plot.subtitle = element_blank(),
                                plot.caption = element_blank())

## Make fancy volcano plot and save
volcano_transcript = EnhancedVolcano(transcript,
                                x = 'log2FoldChange',
                                y = 'padj',
                                lab = transcript$transcript_name,
                                selectLab = c('TGOLN2-201','DYRK1B-202', 'CLU-211', 'NDUFA10-224'),
                                pCutoff = 0.05,
                                ylab = bquote(~-Log[10] ~ italic(Q)),
                                FCcutoff = 1,
                                pointSize = 3.0,
                                colAlpha = 0.6,
                                drawConnectors = TRUE,
                                labSize = 4.5,
                                labCol = 'black',
                                labFace = 'bold',
                                widthConnectors = 0.5,
                                boxedLabels = TRUE,
                                lengthConnectors = unit(0.01, 'npc'),
                                colConnectors = 'black'
                                ) +
                                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                plot.title = element_blank(),
                                plot.subtitle = element_blank(),
                                plot.caption = element_blank())

print(volcano_gene)
print(volcano_transcript)

dev.off()
