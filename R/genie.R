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

library("org.Hs.eg.db")
library(data.table)
library(GENIE3)

data=read.table("../outputs/log-counts.tsv", row.names=1,	header=TRUE, sep="\t", stringsAsFactors=FALSE)

ensembl <- row.names(data)
ensembl <- sub("\\.\\d+$","",ensembl)
row.names(data) <- ensembl

genes <- c("ENSG00000262543","ENSG00000174279","ENSG00000177551","ENSG00000164099","ENSG00000106536",
           "ENSG00000128710","ENSG00000181085","ENSG00000254656","ENSG00000139292","ENSG00000235436",
           "ENSG00000137868","ENSG00000258949","ENSG00000223403","ENSG00000238113","ENSG00000148795",
           "ENSG00000105357","ENSG00000198739","ENSG00000198785","ENSG00000249853","ENSG00000151779",
           "ENSG00000172568","ENSG00000168298","ENSG00000109846","ENSG00000182575","ENSG00000160808",
           "ENSG00000164309","ENSG00000139211","ENSG00000275993","ENSG00000164176","ENSG00000150991",
           "ENSG00000237857","ENSG00000197099","ENSG00000253347","ENSG00000274383","ENSG00000280727",
           "ENSG00000225527","ENSG00000269867","ENSG00000176716","ENSG00000274370","ENSG00000261888")

data <- data.matrix(data[genes,])

row.names(data) <- mapIds(x = org.Hs.eg.db,
                      keys = row.names(data),
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# run one of the next two instructions. Not both!
data <- data[1:nrow(data),1:20] # consider only the SS samples
data <- data[1:nrow(data),21:ncol(data)] # consider only the SL samples

head(data)

regulators <- c("EVX2","NHLH2","POU6F2","HOXD10","LINC01410","MEG9","LOC101927418","DPY19L2P4")
targets <- c("CYP17A1","MYH14","GRIN3A","HS3ST5","NBAS","HIST1H1E","CRYAB","MYL3","CMYA5","AMIGO2","EDIL3",
             "UBC","MAPK15","LGR5","PRSS12","STRA6","RTL1","SMIM28","LRRTM3","FNDC9","SIK1B","NXPH3")

set.seed(123)
res <- GENIE3(data, regulators = regulators, targets = targets)
#write.table(res, file = "../outputs/genie3-1.csv", sep = '\t')

w <- getLinkList(res, threshold = 0.20)
write.table(w, "../outputs/genie3/genie3-ss.tsv", row.names = FALSE, sep = '\t')
