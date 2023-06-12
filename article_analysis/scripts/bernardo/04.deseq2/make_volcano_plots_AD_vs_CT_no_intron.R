#!/usr/bin/env Rscript

## Load libraries
library(ggplot2)
library(EnhancedVolcano)

## Load datasets for differential expression results in genes with and without covariates
gene <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/genes_annotated_results_AD_vs_CT_no_intron.tsv', sep="\t") 

## Load datasets for differential expression results in transcripts with and without covariates
transcript <- read.csv(file = '../../../data/bernardo/processed/04.deseq2/annotated_multiple_transcripts_results_AD_vs_CT.tsv', sep="\t")


# PDF device
mm = (1/25.4)
pdf("../../../figures/bernardo/04.deseq2/volcano_plot_AD_no_intron.pdf", width=(70*mm), height=(90*mm))

## Make fancy volcano plot and save
volcano_gene = EnhancedVolcano(gene,
                                x = 'log2FoldChange',
                                y = 'padj',
                                lab = gene$gene_name,
                                selectLab = c('S100A13','TNFSF12', 'TFG', 'NSD2', 'GABBR1', 'BambuGene290099', 'NEDD9', 'UGP2', 'NAPB'),
                                pCutoff = 0.05,
                                ylab = bquote(~-Log[10] ~ italic(Q)),
                                FCcutoff = 1,
                                pointSize = 0.5,
                                colAlpha = 0.6,
                                drawConnectors = TRUE,
                                axisLabSize = 6,
                                captionLabSize = 6,
                                subtitleLabSize = 6,
                                legendLabSize = 6,
                                legendIconSize = 1,
                                borderWidth = 0.4,
                                titleLabSize = 6,
                                ylim = c(0,11),
                                labSize = 0.9,
                                labCol = 'black',
                                labFace = 'bold',
                                widthConnectors = 0.4,
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
                                axis.ticks=element_line(linewidth=0.4),
                                text = element_text(size = 6))

## Make fancy volcano plot and save
volcano_transcript = EnhancedVolcano(transcript,
                                x = 'log2FoldChange',
                                y = 'padj',
                                lab = transcript$transcript_name,
                                selectLab = c('S100A13-205','TNFSF12-219', 'TNFSF12-203', 'TFG-212', 'NSD2-228', 'GABBR1-207', 'NEDD9-202', 'UGP2-211', 'NAPB-203', 'NAPB-207', "BambuTx2703"),
                                pCutoff = 0.05,
                                ylab = bquote(~-Log[10] ~ italic(Q)),
                                FCcutoff = 1,
                                pointSize = 0.5,
                                colAlpha = 0.6,
                                drawConnectors = TRUE,
                                axisLabSize = 6,
                                captionLabSize = 6,
                                subtitleLabSize = 6,
                                legendLabSize = 6,
                                legendIconSize = 1,
                                borderWidth = 0.4,
                                titleLabSize = 6,
                                ylim = c(0,11),
                                labSize = 0.9,
                                labCol = 'black',
                                labFace = 'bold',
                                widthConnectors = 0.4,
                                boxedLabels = TRUE,
                                maxoverlapsConnectors = Inf,
                                lengthConnectors = unit(0.01, 'npc'),
                                colConnectors = 'black') +
                                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                plot.title = element_blank(),
                                plot.subtitle = element_blank(),
                                plot.caption = element_blank(),
                                axis.ticks=element_line(linewidth=0.4),
                                text = element_text(size = 6))

print(volcano_gene)
print(volcano_transcript)

dev.off()
