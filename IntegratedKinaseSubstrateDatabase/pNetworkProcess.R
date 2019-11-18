source("R/kinase_substrate_databaseIntegration/uniprot.XREF.R")

dat <- readLines("Dat/phosphoNetwork/highResolutionNetwork.csv")
dat <- split(dat, cumsum(grepl("^>", dat)))
seqs <- lapply(dat, function(x) {
  protein <- gsub("^>", "", x[1])
  u <- sapply(USE.NAMES = FALSE, x[-1], function(xx) {
    strsplit(xx, "\t")[[1]]
  })
  u <- t(u)
  aa <- substr(u[, 2], 1, 1)
  position <- as.integer(substr(u[, 2], 2, nchar(u[, 2])))
  r <- cbind(protein, u, aa, position)
  colnames(r) <- c("protein", "flankseq7", "AApos", "kinase", "score", "aminoAcid", "position")
  r
})
seqs <- do.call("rbind", seqs)
seqs[, 2] <- gsub("-", "_", seqs[, 2])


seqs <- cbind(seqs, ACC = gnToAcc(seqs[, "protein"]))
seqs <- cbind(seqs, ID = accTo(sapply(strsplit(seqs[, "ACC"], ";"), "[", 1), to="UniProtKB-ID"))

write.table(seqs, file = "Res/20160803_kinaseSubstrateDatabaseParse/phosphoNetwork.selftab.txt",
            col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")










