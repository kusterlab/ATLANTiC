
source("R/kinase_substrate_databaseIntegration/uniprot.XREF.R")

library(seqinr)

fa <- read.fasta("Dat/Uniprot/uniprot_sprot.fasta", seqtype = "AA")

unip <- read.delim(file = "Res/20160803_kinaseSubstrateDatabaseParse/uniprot.parse.intermid.txt", 
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ids <- paste("sp", unip$proteinAcc, unip$protein, sep = "|")
seq <- sapply(1:length(ids), function(i) {
  seq <- c(rep("_", 15), fa[[ids[i]]], rep("_", 15))
  s <- seq[unip$site[i]:(unip$site[i]+30)]
  c(paste(s, collapse = ""), s[16])
})
seq <- t(seq)
colnames(seq) <- c("flankseq15", "AAseq")
unip <- cbind(unip, seq)
table(unip$aminoAcid == unip$AAseq)
unip <- unip[unip$aminoAcid == unip$AAseq, ]

geneName <- accTo(unip[, "proteinAcc"], to="Gene_Name")
unip <- cbind(unip, geneName)
unip <- data.frame(unip)

write.table(unip, file = "Res/20160803_kinaseSubstrateDatabaseParse/uniprot.selftab.txt",
            col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")





