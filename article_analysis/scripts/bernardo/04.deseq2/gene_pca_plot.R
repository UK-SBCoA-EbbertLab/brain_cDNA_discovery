#!/usr/bin/Rscript

## Import Libraries
library("DESeq2")
library("ggplot2")
#library("svglite")



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

## Filter by only keepting transcripts with XXX or more counts in at least 3 samples
keep_10 <- rowSums(counts(dds) >= 10) >= 3
keep_25 <- rowSums(counts(dds) >= 25) >= 3
keep_50 <- rowSums(counts(dds) >= 50) >= 3
keep_100 <- rowSums(counts(dds) >= 100) >= 3


dds_10 <- dds[keep_10,]
dds_25 <- dds[keep_25,]
dds_50 <- dds[keep_50,]
dds_100 <- dds[keep_100,]


## Create transformed values to make PCA plot
vsd <- vst(dds, blind=FALSE)
vsd_10 <- vst(dds_10, blind=FALSE)
vsd_25 <- vst(dds_25, blind=FALSE)
vsd_50 <- vst(dds_50, blind=FALSE)
vsd_100 <- vst(dds_100, blind=FALSE)


## Calculate things to make PCA plot
pcaData <- plotPCA(vsd, intgroup=c("condition", "sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData_10 <- plotPCA(vsd_10, intgroup=c("condition", "sex"), returnData=TRUE)
percentVar_10 <- round(100 * attr(pcaData_10, "percentVar"))

pcaData_25 <- plotPCA(vsd_25, intgroup=c("condition", "sex"), returnData=TRUE)
percentVar_25 <- round(100 * attr(pcaData_25, "percentVar"))

pcaData_50 <- plotPCA(vsd_50, intgroup=c("condition", "sex"), returnData=TRUE)
percentVar_50 <- round(100 * attr(pcaData_50, "percentVar"))

pcaData_100 <- plotPCA(vsd_100, intgroup=c("condition", "sex"), returnData=TRUE)
percentVar_100 <- round(100 * attr(pcaData_100, "percentVar"))


## Plot PCAs
unfiltered_plot = ggplot(pcaData, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=3) +  xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme(
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill='transparent'),
         legend.box.background = element_rect(fill='transparent')
       )

ten_plot = ggplot(pcaData_10, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=3) +  xlab(paste0("PC1: ",percentVar_10[1],"% variance")) +  ylab(paste0("PC2: ",percentVar_10[2],"% variance")) +
  coord_fixed() + theme(
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill='transparent'),
         legend.box.background = element_rect(fill='transparent')
       )

twenty_five_plot = ggplot(pcaData_25, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=3) +  xlab(paste0("PC1: ",percentVar_25[1],"% variance")) +  ylab(paste0("PC2: ",percentVar_25[2],"% variance")) +
  coord_fixed() + theme(
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill='transparent'),
         legend.box.background = element_rect(fill='transparent')
       )

fifty_plot = ggplot(pcaData_50, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=1) +  xlab(paste0("PC1: ",percentVar_50[1],"% variance")) +  ylab(paste0("PC2: ",percentVar_50[2],"% variance")) +
  coord_fixed() + theme(
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill='transparent'),
         legend.box.background = element_rect(fill='transparent')
       )

hundred_plot = ggplot(pcaData_100, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=1.5) +  xlab(paste0("PC1: ",percentVar_100[1],"% variance")) +  ylab(paste0("PC2: ",percentVar_100[2],"% variance")) +
  coord_fixed() + theme(
         text = element_text(size=5),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         plot.margin = unit(c(0,0,0,0), "mm"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.border = element_rect(color = "black", fill=NA, linewidth=1),
         legend.background = element_rect(fill='transparent', color=NA),
         legend.box.background = element_rect(fill='transparent', color=NA),
         legend.text = element_text(color="black", size=4),
         legend.direction = "horizontal",
         legend.position="top",
         legend.justification = "left",
         legend.spacing.x = unit(0.2, "mm"),
         legend.spacing.y = unit(0.05, "mm"),
         legend.key.height = unit(0, "mm"),
         legend.key.width = unit(4, "mm"),
         legend.margin = margin(-3, 2, -3, 2, "mm"),
         legend.key=element_blank()
       )

ggsave("../../../figures/bernardo/02.quality_control/pca_unfiltered_gene.png", unfiltered_plot, bg="transparent", dpi=1200, units="mm", width=59, height=64)
ggsave("../../../figures/bernardo/02.quality_control/pca_10_gene.png", ten_plot, bg="transparent", dpi=1200, units="mm", width=59, height=64)
ggsave("../../../figures/bernardo/02.quality_control/pca_25_gene.png", twenty_five_plot, bg="transparent", dpi=1200, units="mm", width=59, height=64)
ggsave("../../../figures/bernardo/02.quality_control/pca_50_gene.png", fifty_plot, bg="transparent", dpi=1200, units="mm", width=59, height=64)
ggsave("../../../figures/bernardo/02.quality_control/pca_100_gene.png", hundred_plot, bg="transparent", dpi=1200, units="mm", width=59, height=64)


