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

library(DOSE)
library(enrichplot)
library("AnnotationDbi")
library("org.Hs.eg.db")

# 1) Over-represtation analysis - Up-regulated genes
#---------------------------------------------------

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

x1 <- enrichDO(gene = gene, ont = "DO", pvalueCutoff = 0.06, pAdjustMethod = "BH", minGSSize = 1, maxGSSize = 500,
              qvalueCutoff = 0.2, readable = TRUE)

#x1 <- setReadable(x1, 'org.Hs.eg.db', 'ENTREZID') 
#as.data.frame(x)

# Visualization
p11 <- barplot(x1, showCategory = 20)
#p1 <- cnetplot(x, foldChange = geneList, showCategory = 10)
#p1 <- cnetplot(x, categorySize="pvalue", foldChange = geneList, showCategory = 10)
p12 <- cnetplot(x1, foldChange = geneList, circular=TRUE, colorEdge=TRUE, showCategory = 20)

# 2) Over-represtation analysis - Down-regulated genes
#-----------------------------------------------------

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

x2 <- enrichDO(gene = gene, ont = "DO", pvalueCutoff = 0.06, pAdjustMethod = "BH", minGSSize = 1, maxGSSize = 500,
              qvalueCutoff = 0.2, readable = TRUE)

#x2 <- setReadable(x2, 'org.Hs.eg.db', 'ENTREZID') 

# Visualization
p21 <- barplot(x2, showCategory = 20)
p22 <- cowplot::cnetplot(x2, foldChange = geneList, circular=TRUE, colorEdge=TRUE, showCategory = 20, node_label="gene")

cowplot::plot_grid(p11, p21, ncol=2, labels=LETTERS[1:2], rel_widths=c(1.2, 1.2))
cowplot::plot_grid(p12, p22, ncol=2, labels=LETTERS[1:2], rel_widths=c(1.4, 1.4))

# 3) Perform GSEA. Need the whole dataset(Up + down reg)
#-------------------------------------------------------

data <- read.table("../outputs/0.05.csv", header = TRUE, row.names = 1, sep=",")
data$entrez <- mapIds(x = org.Hs.eg.db,
                       keys = row.names(data),
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
geneList <- data[,2]
names(geneList) <- as.integer(data$entrez)
geneList <- sort(geneList, decreasing = TRUE)
gene <- names(geneList)[abs(geneList) > 1.5]

geneList <- geneList[!is.na(names(geneList))] # remove columns with name is NA
x <- gseDO(geneList, nPerm = 10000, pvalueCutoff = 1.0, minGSSize = 1, maxGSSize = 10)
dim(x)
as.data.frame(x)

# Visualization
enrichplot::dotplot(x)
emapplot(x)
