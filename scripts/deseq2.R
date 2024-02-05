#Install packages if needed:
#install.packages("htmltools")
#BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
library(htmltools)

# Customize by run. Input files are comparison data and comparison metadata
# See the DESeq2 tutorial for information on input files:
# https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# This script is based upon the scripts provided in the above tutorial.

# User must customize:
plotPath = "your/output/folder/path"
compName = "yourComparisonName"
inputPath = "your/input/file/path"

# Make the out file directory
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

# Import comparison data
# File encoding prevents column header issues
countData <- read.csv(paste0(inputPath, compName, '.csv'), header = TRUE, sep = ",", fileEncoding="UTF-8-BOM")
metaData <- read.csv(paste0(inputPath, compName, '_metadata.csv'), header = TRUE, sep = ",", fileEncoding="UTF-8-BOM")

# Convert metadata design column to factor
metaData$genotype = as.factor(metaData$genotype)

# Look at both
head(countData)
head(metaData)

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~genotype, tidy = TRUE)

# Look at the data
dds
dds <- DESeq(dds)

res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table
summary(res) #summary of results
res <- res[order(res$padj),]
head(res)

# Export the results
write.csv(res, paste0(plotPath, "/", "DEresults_", compName,".csv"))


#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))

# User must customize:
metadatacolumn = "user must insert a column name"
gene1 = "user must pick a gene id"
gene2 = "user must pick a gene id"
gene3 = "user must pick a gene id"
gene4 = "user must pick a gene id"
gene5 = "user must pick a gene id"
gene6 = "user must pick a gene id"
plotCounts(dds, gene=gene1, intgroup=metadataColumn)
plotCounts(dds, gene=gene2, intgroup=metadataColumn)
plotCounts(dds, gene=gene3, intgroup=metadataColumn)
plotCounts(dds, gene=gene4, intgroup=metadataColumn)
plotCounts(dds, gene=gene5, intgroup=metadataColumn)
plotCounts(dds, gene=gene6, intgroup=metadataColumn)


# Reset par
par(mfrow=c(1,1))

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3), ylim=c(0,50)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#Transform the raw count data
#vst function will perform variance stabilizing transformation
vsdata <- vst(dds, blind=FALSE)

png(
  paste0(plotPath, "/", "PCAplot_", compName, ".jpg"),
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2
)
plotPCA(vsdata, intgroup=metadatacolumn)
dev.off()
