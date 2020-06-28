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

library("AnnotationDbi")
library("org.Hs.eg.db")
library("magrittr")
library("clusterProfiler")

# 1) Over representation analysis - Upregulated genes
#----------------------------------------------------

topUpDEGs <- read.table("../outputs/top_up.csv", header = TRUE, row.names = 1, sep=",")
topUpDEGs$entrez <- mapIds(x = org.Hs.eg.db,
                       keys = row.names(topUpDEGs),
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
geneList <- topUpDEGs[,2]
names(geneList) <- as.integer(topUpDEGs$entrez)
geneList <- sort(geneList, decreasing = TRUE)
gene <- names(geneList)[abs(geneList) > 1.5]

# GMT processing
wpgmtfile <- system.file("extdata/wikipathways-20200510-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)
wpid2name <- wp2gene %>% dplyr::select(wpid, name)

# Hypergeometric test
ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
as.data.frame(ewp)

# 2) Over representation analysis - Downregulated genes
#------------------------------------------------------

topDw <- read.table("../outputs/top_down.csv", header = TRUE, row.names = 1, sep=",")
topDw$entrez <- mapIds(x = org.Hs.eg.db,
                       keys = row.names(topDw),
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
geneList <- topDw[,2]
names(geneList) <- as.integer(topDw$entrez)
geneList <- sort(geneList, decreasing = TRUE)
gene <- names(geneList)[abs(geneList) > 1.5]

# GMT processing
wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)
wpid2name <- wp2gene %>% dplyr::select(wpid, name)

# Hypergeometric test
ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
as.data.frame(ewp)

# 3) Enrichment analysis
#-----------------------

ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose = FALSE)
as.data.frame(ewp2)


# 4) MSigDb analysis
#-------------------

m_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em <- enricher(gene, TERM2GENE = m_t2g)
em2 <- GSEA(geneList, TERM2GENE = m_t2g)
head(em)
