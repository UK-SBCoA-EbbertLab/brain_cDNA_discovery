#!/usr/bin/Rscript

## Import Libraries
library("DESeq2")
library("ggplot2")
library("dplyr")
library("gridExtra")
library("grid")


## Customize plotPCA function to return more PCs
plotPCA_more_pcs = function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])

  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }


  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], PC4=pca$x[,4], PC5=pca$x[,5], PC6=pca$x[,6], PC7=pca$x[,7], PC8=pca$x[,8], PC9=pca$x[,9], PC10=pca$x[,10], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:10]
    return(d)
  }

  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
        coord_fixed()
}


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
pcaData <- plotPCA_more_pcs(vsd, intgroup=c("condition", "sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
PCs = c(1:10)

df_scree = data.frame(PCs, percentVar)

pdf("../../../figures/bernardo/02.quality_control/pca_scree_plot_gene.pdf", bg="transparent", width=(2.252*2.1), height=2.559)

## Make scree plot
ggplot(data=df_scree, aes(x=PCs, y=percentVar)) + 
  geom_line() + geom_point() +
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ylim(0, 40) + 
  scale_x_continuous(breaks = pretty(df_scree$PCs, n = 10)) +
  theme(text = element_text(size=6))
dev.off()

pdf("../../../figures/bernardo/02.quality_control/PCA_plots_all_pcs_gene.pdf", bg="transparent", width=2.252, height=2.559)

## Plot PC1 vs PC2
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=sex)) +  geom_point(size=1.5) +  xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + coord_fixed(ratio=1.35) + theme(
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
         ) + ylim(-35, 45) + xlim(-50, 60)

## Plot PC1 vs PC3
ggplot(pcaData, aes(PC1, PC3, color=condition, shape=sex)) +  geom_point(size=1.5) +  xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
  coord_fixed() + coord_fixed(ratio=1.35) + theme(
         text = element_text(size=6),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         plot.margin = grid::unit(c(-10, 1,-10,1), "mm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(color = "black", fill=NA, linewidth=1),
         legend.position = "none"
         ) + ylim(-35, 45) + xlim(-50, 60)


## Plot PC1 vs PC4
ggplot(pcaData, aes(PC1, PC4, color=condition, shape=sex)) +  geom_point(size=1.5) +  xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC4: ",percentVar[4],"% variance")) +
  coord_fixed() + coord_fixed(ratio=1.35) + theme(
         text = element_text(size=6),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         plot.margin = grid::unit(c(-10, 1,-10,1), "mm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(color = "black", fill=NA, linewidth=1),
         legend.position = "none"
         ) + ylim(-35, 45) + xlim(-50, 60)


## Plot PC2 vs PC3
ggplot(pcaData, aes(PC2, PC3, color=condition, shape=sex)) +  geom_point(size=1.5) +  xlab(paste0("PC2: ",percentVar[2],"% variance")) +  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
  coord_fixed() + coord_fixed(ratio=1.35) + theme(
         text = element_text(size=6),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         plot.margin = grid::unit(c(-10, 1,-10,1), "mm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(color = "black", fill=NA, linewidth=1),
         legend.position = "none"
         ) + ylim(-35, 45) + xlim(-50, 60)


## Plot PC2 vs PC4
ggplot(pcaData, aes(PC2, PC4, color=condition, shape=sex)) +  geom_point(size=1.5) +  xlab(paste0("PC2: ",percentVar[2],"% variance")) +  ylab(paste0("PC4: ",percentVar[4],"% variance")) +
  coord_fixed() + coord_fixed(ratio=1.35) + theme(
         text = element_text(size=6),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         plot.margin = grid::unit(c(-10, 1,-10,1), "mm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(color = "black", fill=NA, linewidth=1),
         legend.position = "none"
         ) + ylim(-35, 45) + xlim(-50, 60)

## Plot PC3 vs PC4
ggplot(pcaData, aes(PC3, PC4, color=condition, shape=sex)) +  geom_point(size=1.5) +  xlab(paste0("PC3: ",percentVar[3],"% variance")) +  ylab(paste0("PC4: ",percentVar[4],"% variance")) +
  coord_fixed() + coord_fixed(ratio=1.35) + theme(
         text = element_text(size=6),
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         plot.margin = grid::unit(c(-10, 1,-10,1), "mm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(color = "black", fill=NA, linewidth=1),
         legend.position = "none"
       ) + ylim(-35, 45) + xlim(-50, 60)

dev.off()

