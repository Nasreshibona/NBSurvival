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
library("org.Hs.eg.db")

# search_kegg_organism('hsa', by='kegg_code')

# 1) The following is for top upregulated genes. Try for down regulated genes (../outputs/to_down.csv)

data <- read.table("../outputs/top_up.csv", header = TRUE, row.names = 1, sep=",")
data$entrez <- mapIds(x = org.Hs.eg.db,
                      keys = row.names(data),
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")
geneList <- data[,2]
names(geneList) <- as.integer(data$entrez)
geneList <- geneList[!is.na(names(geneList))] # remove columns with name is NA
geneList <- sort(geneList, decreasing = TRUE)
gene <- names(geneList)[abs(geneList) > 1.5]

kk <- enrichKEGG(gene = gene, organism = 'hsa', pvalueCutoff = 0.5)

# 2) GENE SET ENRICHMENT (GSE). HERE WE SHOULD USE THE WHOLE DATASET (UP + DOWN REG)

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
kk2 <- gseKEGG(geneList = geneList, organism = 'hsa', nPerm = 1000, minGSSize = 1, pvalueCutoff = 1.0, verbose = FALSE)
