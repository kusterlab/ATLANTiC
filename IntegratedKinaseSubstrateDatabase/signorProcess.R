source("R/kinase_substrate_databaseIntegration/uniprot.XREF.R")
library(seqinr)

dat <- read.delim("Dat/signor/all_data.tsv")
dat <- dat[dat$EffectMechanism == "phosphorylation" & 
             dat$MechanismResidues != "" &
             dat$MechanismSequences != "", ]
colnames(dat)
dat <- dat[, c("EntityA", "TypeA", "IdA", "EntityB", "TypeB", "IdB", "MechanismResidues", "MechanismSequences")]
dat <- unique(dat)
dat[1:ncol(dat)] <- lapply(dat, as.character)


## manually correct mistakes
dat$MechanismSequences[dat$MechanismSequences == "ERALTEDsTQTSDTA;RALTEDStQTSDTAT;LTEDSTQtSDTATNS;DSTQTSDtATNSTLP;TQTSDTAtNSTLPSA;SDTATNStLPSAEVE"] <- 
  "ERALTEDsTQTSDTA;RALTEDStQTSDTAT;LTEDSTQtSDTATNS;TEDSTQTsDTATNST;DSTQTSDtATNSTLP;TQTSDTAtNSTLPSA;SDTATNStLPSAEVE"
dat$MechanismResidues[ dat$MechanismSequences == "PNEPVSDyINANIIM;NSKPKKSyIATQGCL" ] <- "Tyr304;Tyr327"



dat2 <- lapply(1:nrow(dat), function(i) {
  x <- dat[i, ]
  MechanismResidues <- strsplit(x$MechanismResidues, ";|:")[[1]]
  lr <- length(MechanismResidues)
  if (lr == 1)
    return(x)
  MechanismSequences <- strsplit(x$MechanismSequences, ";")[[1]]
  
  if (lr == length(MechanismSequences)) {
    x <- cbind(matrix(rep(x[1:6], lr), nrow = lr, byrow = TRUE),
               MechanismResidues, MechanismSequences)
    print(dim(x))
  } else 
    x <- x
  colnames(x) <- colnames(dat)
  x
})

dat2 <- do.call(rbind, c(dat2, make.row.names = FALSE))
dat2[sapply(dat2, is.list)] <- lapply(dat2[sapply(dat2, is.list)], unlist)

aa <- substr(dat2$MechanismResidues, 1, 3)
aa[aa == "Thr"] <- "T"
aa[aa == "Ser"] <- "S"
aa[aa == "Tyr"] <- "Y"
pos <- as.integer(substr(dat2$MechanismResidues, 4, nchar(dat2$MechanismResidues)))

dat2$aminoAcid <- aa
dat2$position <- pos
dat2$ID <- accTo(dat2$IdB, to = "UniProtKB-ID")

write.table(dat2, file = "Res/20160803_kinaseSubstrateDatabaseParse/signor.selftab.txt",
            col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")

