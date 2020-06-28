library("org.Hs.eg.db")
library(data.table)
library(GENIE3)
library(edgeR)
library(DESeq2)

# Raw-counts
data=read.table("../outputs/raw-counts.csv", row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

ensembl <- row.names(data)
ensembl <- sub("\\.\\d+$","",ensembl)
row.names(data) <- ensembl

group <- factor(c(rep("S",20), rep("L",12)))
dge <- DGEList(counts = data, group = group)
keep <- filterByExpr(dge)
summary(keep)
dge <- dge[keep, , keep.lib.sizes=FALSE]
genes <- rownames(dge$counts)

# Log-counts
data = read.table("../outputs/log-counts.tsv",row.names=1,header=TRUE,sep="\t",stringsAsFactors=FALSE)
ensembl <- row.names(data)
ensembl <- sub("\\.\\d+$","",ensembl)
row.names(data) <- ensembl
data <- data[genes, ]

data <- data[1:nrow(data),1:20]
data <- data.matrix(data)

dim(data)
head(data, 2)

set.seed(123)
res <- GENIE3(data, regulators = NULL, targets = NULL, nCores=1)
w <- getLinkList(res)
write.table(w, "../genie3-ss.tsv", row.names = FALSE, sep = '\t')

