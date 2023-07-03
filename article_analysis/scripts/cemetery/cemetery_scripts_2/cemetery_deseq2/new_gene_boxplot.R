#!/usr/bin/Rscript

## Import Libraries
library("DESeq2")
library("ggplot2")



## Import counts matrix
cts <- as.matrix(read.csv("../../../data/bernardo/processed/04.deseq2/transcript_counts_unfiltered.tsv", sep="\t", row.names="transcript_id"))

## Import metadata
coldata <- read.csv("../../../data/bernardo/processed/04.deseq2/experimental_design.tsv", sep="\t", row.names=1)

coldata$condition <- factor(coldata$condition)
coldata$sex <- factor(coldata$sex)


## Make sure order of columns in metadata and counts matrix is the same... TRUE = good, FALSE = bad
all(rownames(coldata) == colnames(cts))

any("BambuTx2073" == rownames(cts))


## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)


# PDF device
pdf("../../../figures/bernardo/03.gene_and_transcripts_descriptive_stats/graduate_showcase_poster/THE_new_transcript_expression.pdf")

for (transcript_id in c("BambuTx2703")) {

    print(transcript_id)
    d = plotCounts(dds, gene=transcript_id, intgroup="condition", returnData=TRUE, normalized=TRUE) 

    plot_title = "GL000214.1: 125344-133145 (-)"

    plot = ggplot(d, aes(x=condition, y=count, fill=condition)) + ylab("DESeq2 Normalized Counts") +
    geom_boxplot(alpha=0.5, outlier.shape=NA)+ ggtitle(plot_title) + theme(text = element_text(size = 24), plot.title = element_text(color="red", size=24, face="bold.italic", hjust = 0.5), legend.position="none") +
    geom_point(position=position_jitter(w=0.1,h=0))
    
    print(plot)
}

dev.off()


