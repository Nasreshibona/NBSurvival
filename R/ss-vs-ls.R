library(edgeR)
library(DESeq2)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# raw-counts.csv created using the script ../python/raw-counts.ipynb
# change directory to the R folder
data=read.table("../outputs/raw-counts.csv", row.names=1,	header=TRUE, sep="\t", stringsAsFactors=FALSE)

group <- factor(c(rep("S",20), rep("L",12)))
dge <- DGEList(counts = data, group = group)
keep <- filterByExpr(dge)
summary(keep)
dge <- dge[keep, , keep.lib.sizes=FALSE]
m <- dge$counts

pData <- read.table("../outputs/pData.txt", row.names=1, header=TRUE, sep="\t")
all(rownames(pData)==colnames(m))
metadata <- data.frame(labelDescription=c("Condition Name", "Sample ID"), row.names = c("conditionName", "sampleId"))
phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata=metadata)
conditionName <- phenoData$conditionName
eset <- ExpressionSet(assayData = m, phenoData = phenoData)
colData <- DataFrame(pData(eset))
se <- SummarizedExperiment(exprs(eset), colData=colData)

mode(assay(se)) <- "integer"

dds <- DESeqDataSet(se, design = ~ conditionName)
dds$conditionName <- relevel(dds$conditionName, ref="L")
dds <- DESeq(dds)

res <- results(dds, lfcThreshold = 0.3, alpha = 0.05, altHypothesis = "greaterAbs")
summary(res)
resultsNames(dds)
mcols(res, use.names = T)

res <- res[order(res$padj),]
resSig <- res[res$padj < 0.05 & !is.na(res$padj), ]
summary(resSig)

# ANNOTATE GENE SYMBOLS
#----------------------

resSig$transcriptId <- rownames(resSig)
ensembl <- row.names(resSig)
ensembl <- sub("\\.\\d+$","",ensembl)
row.names(resSig) <- ensembl
head(resSig)
resSig$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = row.names(resSig),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

write.table(resSig, "../outputs/0.05.csv", sep=",")

# downregulated genes
top_down <- resSig[order(resSig$log2FoldChange, -resSig$baseMean), ][1:19,]
write.csv(top_down, "../outputs/top_down.csv")
# upregulated genes
top_up <- resSig[order(-resSig$log2FoldChange, -resSig$baseMean), ][1:21,]
write.csv(top_up, "../outputs/top_up.csv")


# ***** PLOTTING *********

# Run again this to reset the index
res <- results(dds, lfcThreshold = 0.3, alpha = 0.05, altHypothesis = "greaterAbs")
res <- res[order(res$padj),]
resSig <- res[res$padj < 0.05 & !is.na(res$padj), ]

dds_rlog <- vst(dds, blind = FALSE)
plotPCA(dds_rlog, intgroup = "conditionName", ntop = 2000) +
  theme_bw() + # remove default ggplot2 theme
  geom_point(size = 5) + # Increase point size
  scale_y_continuous(limits = c(-50, 50)) + # change limits to fix figure dimensions
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 500 most variable genes")

mat <- assay(dds_rlog[row.names(resSig)])[1:20, ]

annotation_col = data.frame(
  Group = factor(colData(dds_rlog)$conditionName),
  row.names = colData(dds_rlog)$sampleid
)

rownames(annotation_col) <- colnames(mat)

pheatmap(mat = mat, 
         scale = "row", # Scale genes to Z-score (how many standard deviations)
         annotation_col = annotation_col, # Add multiple annotations to the samples
         fontsize = 6.5, # Make fonts smaller
         show_colnames = T)
