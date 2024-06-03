## RNA-seq analysis with DESeq2

# Set working directory
# cd /home/zcqsjwi/Scratch/counts
setwd("~/Documents/Rotation2/CLOZvsCTRLvsHAL")

# Import the data
#countdata <- read.table("counts.txt", header=TRUE, row.names=1)
countdata <- read.delim("~/Documents/Rotation2/countdata.txt", header=TRUE, comment.char="#", row.names=1)

# Only need the counts so trim tables 1st 5 columns
countdata <- countdata[-c(1:5)]

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

# Remove .......aligned_reads. from column names
colnames(countdata) <- gsub(".........aligned_reads.", "", colnames(countdata))

# Remove anything after underscore in column names
colnames(countdata) <- gsub("_.*", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition (first four are controls, second four and third four contain two different experiments)
(condition <- factor(c(rep("cloz", 10), rep("ctrl", 9), rep("hal", 10))))

# Analysis with DESeq2

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal component analysis
## Could do with built-in DESeq2 function:
DESeq2::plotPCA(rld, intgroup="condition")
## I (Stephen Turner) like mine better:
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    } else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  # rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  # pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  # terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()


# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
## Could do with built-in DESeq2 function:
#DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=1, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-4, 6.5))
dev.off()

## Annotation with gene symbols and filtering
#install.packages("dplyr")
library(dplyr)
library(BiocManager)
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

#Create new data frame with annotation information
anno <- AnnotationDbi::select(org.Mm.eg.db,keys=resdata$Gene,
                              columns=c("SYMBOL","GENENAME"),
                              keytype="ENSEMBL")

#Identify duplicated entries
dup_ids <- anno$ENSEMBL[duplicated(anno$ENSEMBL)]
filter(anno, ENSEMBL %in% dup_ids) %>% 
  arrange(ENSEMBL)

#Select the columns we want
anno <- AnnotationDbi::select(org.Mm.eg.db,keys=resdata$Gene,
                              columns=c("ENSEMBL","SYMBOL","GENENAME","ENTREZID"),
                              keytype="ENSEMBL") %>% 
  filter(!duplicated(ENSEMBL))

#Rename first column in resdata file from Gene to Ensembl
resdata.labelled <- resdata
colnames(resdata.labelled)
names(resdata.labelled)[names(resdata.labelled) == "Gene"] <- "ENSEMBL"

#Bind the annotation information to the results data frame
results.annotated <- left_join(resdata.labelled, anno,by="ENSEMBL")

#Change the order of the columns (move last three columns to start)
results.annotated <- results.annotated[,c(37,38,1,39,2:36)]

#Write the unfiltered annotated results to excel
write.csv(results.annotated, file="diffexpr-results-annotated.csv")

#Filter for abs(Log2FoldChange)>=1 and padj<0.05
results.p0.05 <- filter(results.annotated, padj < 0.05)
results.p0.05.fc1 <- filter(results.p0.05, abs(log2FoldChange) >= 1)
write.csv(results.p0.05.fc1, file="diffexpr-results-annotated-p0.05-fc1.csv")