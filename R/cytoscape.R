library("org.Hs.eg.db")

# Regulatory gene interactions for SS samples. Change filename to filt-genie3-sa.tsv for LS interactions

data=read.table("../outputs/genie3/fromcluster/filt-genie3-sa.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dim(data)
data <- data[data$weight>0.00251, ]
dim(data)

reg <- data$regulatoryGene
regS <- mapIds(x=org.Hs.eg.db, keys=reg, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
sum(is.na(regS))
for (i in 1:length(regS)) {
  if (is.na(regS[i])) {
    regS[i] <- names(regS)[i]
  }
}
sum(is.na(regS))
data$regulatoryGene <- regS

tar <- data$targetGene
tarS <- mapIds(x=org.Hs.eg.db, keys=tar, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
sum(is.na(regS))
for(i in 1:length(tarS)) {
  if (is.na(tarS[i])) {
    tarS[i] <- names(tarS)[i]
  }
}
sum(is.na(tarS))
data$targetGene <- tarS

# Change filename to genie3-ss.tsv to store for LS samples.
write.table(data, "../outputs/genie3/fromcluster/genie3-ls.tsv", row.names = FALSE, sep = '\t')

