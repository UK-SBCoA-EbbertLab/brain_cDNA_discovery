#!/usr/bin/Rscript

## Import Libraries
library("DESeq2")
library("ggplot2")



## Import counts matrix
cts <- as.matrix(read.csv("../../../data/bernardo/processed/04.deseq2/gene_counts_unfiltered.tsv", sep="\t", row.names="gene_id"))

## Import metadata
coldata <- read.csv("../../../data/bernardo/processed/04.deseq2/experimental_design.tsv", sep="\t", row.names=1)

coldata$condition <- factor(coldata$condition)
coldata$sex <- factor(coldata$sex)


## Import gene converter
converter = read.csv("../../../data/bernardo/processed/04.deseq2/AD_gene_converter.tsv", sep="\t")


## Make sure order of columns in metadata and counts matrix is the same... TRUE = good, FALSE = bad
all(rownames(coldata) == colnames(cts))



## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ sex)



# PDF device
pdf("../../../figures/bernardo/04.deseq2/AD_gene_level_plots_M_vs_F.pdf")

for (gene_id in converter$gene_id) {

    d = plotCounts(dds, gene=gene_id, intgroup="sex", returnData=TRUE, normalized=FALSE)

    gene_name = converter[which(converter$gene_id == gene_id), ]$gene_name

    plot = ggplot(d, aes(x=sex, y=count, fill=sex)) +
    geom_boxplot(alpha=0.5, outlier.shape=NA)+ ggtitle(gene_name) + theme(plot.title = element_text(color="red", size=24, face="bold.italic", hjust = 0.5), legend.position="none") +
    geom_point(position=position_jitter(w=0.1,h=0))
    
    print(plot)
}

dev.off()


