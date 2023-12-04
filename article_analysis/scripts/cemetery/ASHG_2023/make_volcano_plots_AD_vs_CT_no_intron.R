#!/usr/bin/env Rscript

## Load libraries
library(ggplot2)
library(EnhancedVolcano)

## Load datasets for differential expression results in genes with and without covariates
gene <- read.csv(file = '../../data/processed/deseq2/genes_annotated_results_AD_vs_CT_no_intron.tsv', sep="\t") 

## Load datasets for differential expression results in transcripts with and without covariates
transcript <- read.csv(file = '../../data/processed/deseq2/annotated_multiple_transcripts_results_AD_vs_CT.tsv', sep="\t")


# PDF device
pdf("../../figures/ASHG_2023/volcano_plot_AD_no_intron.pdf", width=10.5, height=6.2)

## Make fancy volcano plot and save
volcano_gene = EnhancedVolcano(gene,
                                x = 'log2FoldChange',
                                y = 'padj',
                                lab = gene$gene_name,
                                selectLab = c('TNFSF12', 'NSD2', 'UGP2'),
                                pCutoff = 0.05,
                                ylab = bquote(~-Log[10] ~ italic(Q)),
                                FCcutoff = 1,
                                pointSize = 8,
                                colAlpha = 0.6,
                                drawConnectors = TRUE,
                                axisLabSize = 12,
                                captionLabSize = 12,
                                subtitleLabSize = 12,
                                legendLabSize = 12,
                                legendIconSize = 1,
                                borderWidth = 2,
                                titleLabSize = 12,
                                ylim = c(0,11),
                                labSize = 12,
                                labCol = 'black',
                                labFace = 'bold',
                                widthConnectors = 1,
                                boxedLabels = TRUE,
                                maxoverlapsConnectors = Inf,
                                lengthConnectors = unit(0.01, 'npc'),
                                colConnectors = 'black'
                                ) +     
                                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                plot.title = element_blank(),
                                plot.subtitle = element_blank(),
                                plot.caption = element_blank(),
                                axis.ticks=element_line(linewidth=2),
                                text = element_text(size = 12))

## Make fancy volcano plot and save
volcano_transcript = EnhancedVolcano(transcript,
                                x = 'log2FoldChange',
                                y = 'padj',
                                lab = transcript$transcript_name,
                                selectLab = c('TNFSF12-219', 'TNFSF12-203', 'NSD2-228', 'UGP2-211'),
                                pCutoff = 0.05,
                                ylab = bquote(~-Log[10] ~ italic(Q)),
                                FCcutoff = 1,
                                pointSize = 8,
                                colAlpha = 0.6,
                                drawConnectors = TRUE,
                                axisLabSize = 12,
                                captionLabSize = 12,
                                subtitleLabSize = 12,
                                legendLabSize = 12,
                                legendIconSize = 1,
                                borderWidth = 2,
                                titleLabSize = 12,
                                ylim = c(0,11),
                                labSize = 12,
                                labCol = 'black',
                                labFace = 'bold',
                                widthConnectors = 1,
                                boxedLabels = TRUE,
                                maxoverlapsConnectors = Inf,
                                lengthConnectors = unit(0.01, 'npc'),
                                colConnectors = 'black') +
                                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                plot.title = element_blank(),
                                plot.subtitle = element_blank(),
                                plot.caption = element_blank(),
                                axis.ticks=element_line(linewidth=2),
                                text = element_text(size = 12))

print(volcano_gene)
print(volcano_transcript)

dev.off()
