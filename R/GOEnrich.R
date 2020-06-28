# Copyright (C) 2020  Hocine Bendou <hocine@sanbi.ac.za>
#                     Abdulazeez Giwa <3901476@myuwc.ac.za>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>

library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)

# 1) Over representation analysis - Up-regulated genes
#-----------------------------------------------------

topUpDEGs <- read.table("../outputs/dge/top_up.csv", header = TRUE, row.names = 1, sep=",")
topUpDEGs$entrez <- mapIds(x = org.Hs.eg.db,
                           keys = row.names(topUpDEGs),
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multiVals = "first")

geneList <- topUpDEGs[,2]
names(geneList) <- as.integer(topUpDEGs$entrez)
geneList <- geneList[!is.na(names(geneList))] # remove columns with name is NA
geneList <- sort(geneList, decreasing = TRUE)
gene <- names(geneList)[abs(geneList) > 1.5]
gene.df <- bitr(gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)

x11 <- groupGO(gene = gene, OrgDb = org.Hs.eg.db, ont = "MF", level = 3, readable = TRUE)
x11 <- x11[x11$Count>0, ] # Remove terms with count = 0
head(x11)
x12 <- groupGO(gene = gene, OrgDb = org.Hs.eg.db, ont = "CC", level = 3, readable = TRUE)
x12 <- x12[x12$Count>0, ] # Remove terms with count = 0
head(x12)
x13 <- groupGO(gene = gene, OrgDb = org.Hs.eg.db, ont = "BP", level = 3, readable = TRUE)
x13 <- x13[x13$Count>0, ] # Remove terms with count = 0
head(x13)

#x13 <- x13[x13$Count>0 & x13$geneID %like% "PRSS12", ] # TRY IT TO SEE RESULTS FOR ONE GENE

# as.data.frame(x11)

# Enrich GO
x1 <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, readable = TRUE, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
x2 <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, readable = TRUE, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
x3 <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, readable = TRUE, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.08, qvalueCutoff = 0.07)


# FILTER RESULT
x1 <- x1[x1$geneID %like% "EDIL3", ] # TRY IT TO SEE RESULTS FOR ONE GENE

p1 <- barplot(x1, showCategory = 20)
p2 <- barplot(x2, showCategory = 20)
p3 <- barplot(x3, showCategory = 20)
cowplot::plot_grid(p1, p3, ncol=2, labels=LETTERS[1:2], rel_widths=c(1.5, 1.5))

# 2) Over representation analysis - Down-regulated genes
#-------------------------------------------------------

topDwDEGs <- read.table("../outputs/dge/top_down.csv", header = TRUE, row.names = 1, sep=",")
topDwDEGs$entrez <- mapIds(x = org.Hs.eg.db,
                          keys = row.names(topDwDEGs),
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")

geneList <- topDwDEGs[,2]
names(geneList) <- as.integer(topDwDEGs$entrez)
geneList <- geneList[!is.na(names(geneList))] # remove columns with name is NA
geneList <- sort(geneList, decreasing = TRUE)
gene <- names(geneList)[abs(geneList) > 1.5]
gene.df <- bitr(gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)
head(gene.df)


x1 <- groupGO(gene = gene, OrgDb = org.Hs.eg.db, ont = "MF", level = 3, readable = TRUE)
x2 <- groupGO(gene = gene, OrgDb = org.Hs.eg.db, ont = "CC", level = 3, readable = TRUE)
x3 <- groupGO(gene = gene, OrgDb = org.Hs.eg.db, ont = "BP", level = 3, readable = TRUE)

# FILTER RESULT
x1 <- x1[x1$Count>0, ] # Remove terms with count = 0
x2 <- x2[x2$Count>0, ] # Remove terms with count = 0
x3 <- x3[x3$Count>0, ] # Remove terms with count = 0
#x3 <- x3[x3$Count>0 & x3$geneID %like% "HOXD10", ] # TRY IT TO SEE RESULTS FOR ONE GENE
#dim(x3)
#as.data.frame(x3)

# Enrich GO
x1 <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, readable = TRUE, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
x2 <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, readable = TRUE, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.5, qvalueCutoff = 0.5)
x3 <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, readable = TRUE, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.08, qvalueCutoff = 0.07)

# FILTER RESULT
x1 <- x1[x1$geneID %like% "HOXD10", ] # TRY IT TO SEE RESULTS FOR ONE GENE

p1 <- barplot(x1, showCategory = 20)
p2 <- barplot(x2, showCategory = 20)
p3 <- barplot(x3, showCategory = 20)
cowplot::plot_grid(p1, p3, ncol=2, labels=LETTERS[1:2], rel_widths=c(1.5, 1.5))

# 3) GSEA use the whole dataset (UP + Down Reg) 
#----------------------------------------------

data <- read.table("../outputs/0.05.csv", header = TRUE, row.names = 1, sep=",")
data$entrez <- mapIds(x = org.Hs.eg.db,
                      keys = row.names(data),
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")

geneList <- data[,2]
names(geneList) <- as.integer(data$entrez)
geneList <- geneList[!is.na(names(geneList))] # remove columns with name is NA
geneList <- sort(geneList, decreasing = TRUE)

x1 <- gseGO(geneList=geneList, OrgDb=org.Hs.eg.db, ont="MF", nPerm=1000, minGSSize=2, maxGSSize=50, pvalueCutoff=0.5, verbose=FALSE)
x2 <- gseGO(geneList=geneList, OrgDb=org.Hs.eg.db, ont="CC", nPerm=1000, minGSSize=2, maxGSSize=50, pvalueCutoff=0.7, verbose=FALSE)
x3 <- gseGO(geneList=geneList, OrgDb=org.Hs.eg.db, ont="BP", nPerm=500, minGSSize=2, maxGSSize=50, pvalueCutoff=0.5, verbose=FALSE)

x1 <- setReadable(x1, 'org.Hs.eg.db')
x2 <- setReadable(x2, 'org.Hs.eg.db')
x3 <- setReadable(x3, 'org.Hs.eg.db')

p1 <- emapplot(x1)
p2 <- emapplot(x2)
p3 <- emapplot(x3)
