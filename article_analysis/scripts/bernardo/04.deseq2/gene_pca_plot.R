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


## Make sure order of columns in metadata and counts matrix is the same... TRUE = good, FALSE = bad
all(rownames(coldata) == colnames(cts))



## Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

## Create transformed values to make PCA plot
vsd <- vst(dds, blind=FALSE)

## Calculate things to make PCA plot
pcaData <- plotPCA(vsd, intgroup=c("condition", "sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


## Plot PCAs
pca_plot = ggplot(pcaData, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=1.5) +  xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + coord_fixed(ratio=1.65) + theme(
         text = element_text(size=6),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         plot.margin = grid::unit(c(-10, 1,-10,1), "mm"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.border = element_rect(color = "black", fill=NA, linewidth=1),
         legend.background = element_rect(fill='transparent', color="grey83", linewidth=0.5),
         legend.box.background = element_rect(fill='transparent', color=NA),
         legend.text = element_text(color="black", size=5),
         legend.direction = "horizontal",
         legend.position = c(0.975, 0.15),
         legend.justification = "right",
         legend.spacing.x = unit(1, "mm"),
         legend.spacing.y = unit(2, "mm"),
         legend.key.height = unit(2, "mm"),
         legend.key.width = unit(2, "mm"),
         legend.margin = margin(1, 1, 1, 1, "mm"),
         legend.key=element_blank()
       )

ggsave("../../../figures/bernardo/02.quality_control/pca_unfiltered_gene.pdf", pca_plot, bg="transparent", dpi=1200, units="mm", width=57.2, height=63)


