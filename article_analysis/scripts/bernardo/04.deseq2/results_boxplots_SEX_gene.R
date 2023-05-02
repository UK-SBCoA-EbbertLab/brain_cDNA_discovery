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
converter = read.csv("../../../data/bernardo/processed/04.deseq2/gene_results_converter_SEX.tsv", sep="\t")


## Make sure order of columns in metadata and counts matrix is the same... TRUE = good, FALSE = bad
all(rownames(coldata) == colnames(cts))



## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ sex)



# PDF device
pdf("../../../figures/bernardo/04.deseq2/results_SEX_boxplots_gene_level.pdf")

for (gene_id in converter$gene_id) {

    chr = converter$chr[converter$gene_id == gene_id]
    d = plotCounts(dds, gene=gene_id, intgroup="sex", returnData=TRUE, normalized=FALSE)

    gene_name = converter[which(converter$gene_id == gene_id), ]$gene_name

    plot_title = paste(chr, gene_name, sep= " - ")
    plot_title = paste("chr", plot_title, sep=":")

    plot = ggplot(d, aes(x=sex, y=count, fill=sex)) + ylab("DESeq2 normalized counts") + 
    geom_boxplot(alpha=0.5, outlier.shape=NA)+ ggtitle(plot_title) + theme(plot.title = element_text(color="red", size=24, face="bold.italic", hjust = 0.5), legend.position="none") +
    geom_point(position=position_jitter(w=0.1,h=0))
    
    print(plot)
}

dev.off()


