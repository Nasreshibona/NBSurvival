library(data.table)
data=read.table("../outputs/genie3/fromcluster/genie3-ss.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
genes <- c("SMIM28","HOXD10","PRSS12","LGR5","EVX2","NHLH2")
data <- data[data$regulatoryGene %in% genes | data$targetGene %in% genes, ]
data <- data[!startsWith(data$targetGene, "ENSG") & !startsWith(data$regulatoryGene, "ENSG"), ] 
write.table(data, "../outputs/genie3/net-2.tsv", row.names = FALSE, sep = '\t')

# only between SS DEGs Net-1
data=read.table("../outputs/genie3/fromcluster/genie3-ss.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
genes <- c("SMIM28","RNF150","LGR5","GLRA2","HOXD10","TIFAB","PRSS12","PLA2G2D","IRF6","EVX2","ELFN1","NHLH2","LINC01266")
data <- data[data$regulatoryGene %in% genes & data$targetGene %in% genes, ]
write.table(data, "../outputs/genie3/net-so-1.tsv", row.names = FALSE, sep = '\t')

# only between SS DEGs Net-2
data=read.table("../outputs/genie3/fromcluster/genie3-ss.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
genes <- c("MEG9","CRYM","UBC","IDNK","ATP6V1G1","STRA6","GRIN3A","KIF1B","IGFL4","FNDC9","ENSG00000258590","ERN1","NXPH3","MYL3","ECHDC3")
data <- data[data$regulatoryGene %in% genes & data$targetGene %in% genes, ]
write.table(data, "../outputs/genie3/net-so-2.tsv", row.names = FALSE, sep = '\t')

# only between SS DEGs Net-3
data=read.table("../outputs/genie3/fromcluster/genie3-ss.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
genes <- c("HS3ST5","TAGLN3","HIST1H1E")
data <- data[data$regulatoryGene %in% genes & data$targetGene %in% genes, ]
write.table(data, "../outputs/genie3/net-so-3.tsv", row.names = FALSE, sep = '\t')


# only between SS DEGs Net-4
data=read.table("../outputs/genie3/fromcluster/genie3-ss.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
genes <- c("MAPK15","ENSG00000275807","ENSG00000269867","CHRM1","NBAS","BBS5","EDIL3","PVR","CYP17A1")
data <- data[data$regulatoryGene %in% genes & data$targetGene %in% genes, ]
write.table(data, "../outputs/genie3/net-so-4.tsv", row.names = FALSE, sep = '\t')
