#!/usr/bin/Rscript

## Import Libraries
library("DESeq2")
library("ggplot2")



## Import counts matrix
cts <- as.matrix(read.csv("../../../data/bernardo/processed/04.deseq2/new_mito_unique.tsv", sep="\t", row.names="transcript_id"))

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



# PDF device
pdf("../../../figures/bernardo/04.deseq2/new_mito_unique_boxplots.pdf")

for (transcript_id in row.names(cts)) {

    d = plotCounts(dds, gene=transcript_id, intgroup="condition", returnData=TRUE)

    plot_title = transcript_id

    plot = ggplot(d, aes(x=condition, y=count, fill=condition)) +
    geom_boxplot(alpha=0.5, outlier.shape=NA)+ ggtitle(plot_title) + theme(plot.title = element_text(color="red", size=24, face="bold.italic", hjust = 0.5), legend.position="none") +
    geom_point(position=position_jitter(w=0.1,h=0))
    
    print(plot)
}

dev.off()


