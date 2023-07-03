#!/usr/bin/Rscript

## Import Libraries
library("DESeq2")
library("ggplot2")



## Import counts matrix
cts_transcript <- as.matrix(read.csv("../../data/processed/deseq2/transcript_counts_unfiltered.tsv", sep="\t", row.names="transcript_id"))
cts_gene <- as.matrix(read.csv("../../data/processed/deseq2/gene_counts_unfiltered_no_introns.tsv", sep="\t", row.names="gene_id"))


## Import metadata
coldata <- read.csv("../../data/processed/deseq2/experimental_design.tsv", sep="\t", row.names=1)

coldata$condition <- factor(coldata$condition)
coldata$sex <- factor(coldata$sex)


## Import gene converter
converter = read.csv("../../data/processed/deseq2/transcript_results_converter_AD.tsv", sep="\t")


## Make sure order of columns in metadata and counts matrix is the same... TRUE = good, FALSE = bad
all(rownames(coldata) == colnames(cts_transcript))
all(rownames(coldata) == colnames(cts_gene))


## Create DESeq2 object for transcript level data
dds_transcript <- DESeqDataSetFromMatrix(countData = cts_transcript,
                              colData = coldata,
                              design = ~ condition)

## Create DESeq2 object for gene level data
dds_gene <- DESeqDataSetFromMatrix(countData = cts_gene,
                              colData = coldata,
                              design = ~ condition)


## Initiate list of gene ids plotted
list_gene_ids = c()


# PDF device
pdf("../../figures/paper_figures/figure_6/results_AD_boxplots_transcript_level_no_intron.pdf")

## Loop through transcript ids
for (transcript_id in converter$transcript_id) {

    ## Get gene id and chromosome for the current transcript
    gene_id = converter$gene_id[converter$transcript_id == transcript_id]
    chr = converter$chr[converter$transcript_id == transcript_id]

    ## If it is the first transcript for a given gene, then display boxplot for gene level expression
    if (!(gene_id %in% list_gene_ids)) {

        d = plotCounts(dds_gene, gene=gene_id, intgroup="condition", returnData=TRUE, normalized=TRUE)

        gene_name = converter[which(converter$transcript_id == transcript_id), ]$gene_name

        plot_title = paste(chr, gene_name, sep=" - ")
        plot_title = paste("chr", plot_title, sep=":")

        plot = ggplot(d, aes(x=condition, y=count, fill=condition)) + ylab("DESeq2 normalized counts") + 
        geom_boxplot(alpha=0.5, outlier.shape=NA)+ ggtitle(plot_title) + theme(plot.title = element_text(color="red", size=24, face="bold.italic", hjust = 0.5), legend.position="none") +
        geom_point(position=position_jitter(w=0.1,h=0))

        print(plot)

        list_gene_ids = append(list_gene_ids, gene_id)
    }

    ## Plot transcript level data
    d = plotCounts(dds_transcript, gene=transcript_id, intgroup="condition", returnData=TRUE, normalized=TRUE)

    gene_name = converter[which(converter$transcript_id == transcript_id), ]$gene_name

    plot_title = paste(chr, gene_name, transcript_id, sep=" - ") 
    plot_title = paste("chr", plot_title, sep=":")

    plot = ggplot(d, aes(x=condition, y=count, fill=condition)) + ylab("DESeq2 normalized counts") +
    geom_boxplot(alpha=0.5, outlier.shape=NA)+ ggtitle(plot_title) + theme(plot.title = element_text(color="red", size=24, face="bold.italic", hjust = 0.5), legend.position="none") +
    geom_point(position=position_jitter(w=0.1,h=0))

    print(plot)
}

dev.off()

