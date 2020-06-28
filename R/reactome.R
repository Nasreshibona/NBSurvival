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

library(ReactomePA)
library(enrichplot)

# 1) Over-representation (enrichment) analysis - up-regulated genes
#------------------------------------------------------------------

topUpDEGs <- read.table("../outputs/top_up.csv", header = TRUE, row.names = 1, sep=",")
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

x1 <- enrichPathway(gene, organism = 'human', minGSSize = 1, pvalueCutoff = 0.1, qvalueCutoff = 0.1, readable = TRUE)
dim(x1)
as.data.frame(x1)

# Visualization
p1 <- cnetplot(x1, foldChange = geneList, circular=TRUE, colorEdge=TRUE, showCategory = 10)

# 2) Over-representation (enrichment) analysis - down-regulated genes
#--------------------------------------------------------------------

topDwDEGs <- read.table("../outputs/top_down.csv", header = TRUE, row.names = 1, sep=",")
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

x2 <- enrichPathway(gene, organism = 'human', minGSSize = 1, pvalueCutoff = 0.1, qvalueCutoff = 0.1, readable = TRUE)
dim(x2)
as.data.frame(x2)

# Visualization
p2 <- cnetplot(x2, foldChange = geneList, circular=TRUE, colorEdge=TRUE, showCategory = 10)
cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2], rel_widths=c(1.2, 1.2))

# 3) GSEA (use the whole degs i.e., up + down reg genes)
#-------------------------------------------------------

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

gse <- gsePathway(geneList = geneList, organism = "human", minGSSize = 1, pvalueCutoff = 1.0, verbose = FALSE)

gse <- setReadable(gse, 'org.Hs.eg.db')
#gse[!is.na(gse$core_enrichment) & gse$enrichmentScore > 0, ]
gse

# visualization
emapplot(gse)
