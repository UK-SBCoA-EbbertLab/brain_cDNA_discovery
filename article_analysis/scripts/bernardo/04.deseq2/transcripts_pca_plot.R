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


# PDF device
pdf("../../../figures/bernardo/02.quality_control/transcript_level_PCAs.pdf")

## Plot PCAs
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=3) +  xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + ggtitle("Unfiltered")

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=3) +  xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + ggtitle("Unfiltered")

ggplot(pcaData_10, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=3) +  xlab(paste0("PC1: ",percentVar_10[1],"% variance")) +  ylab(paste0("PC2: ",percentVar_10[2],"% variance")) +
  coord_fixed() + ggtitle("10+ counts in 4/12 samples")

ggplot(pcaData_25, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=3) +  xlab(paste0("PC1: ",percentVar_25[1],"% variance")) +  ylab(paste0("PC2: ",percentVar_25[2],"% variance")) +
  coord_fixed() + ggtitle("25+ counts in 4/12 samples")

ggplot(pcaData_50, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=3) +  xlab(paste0("PC1: ",percentVar_50[1],"% variance")) +  ylab(paste0("PC2: ",percentVar_50[2],"% variance")) +
  coord_fixed() + ggtitle("50+ counts in 4/12 samples")

ggplot(pcaData_100, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=3) +  xlab(paste0("PC1: ",percentVar_100[1],"% variance")) +  ylab(paste0("PC2: ",percentVar_100[2],"% variance")) +
  coord_fixed() + ggtitle("100+ counts in 4/12 samples")

dev.off()


