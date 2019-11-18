# xref <- read.delim("Dat/Uniprot/HUMAN_9606_idmapping.dat", header = FALSE, stringsAsFactors=FALSE)
# colnames(xref) <- c("ACC", "source", "value")
# saveRDS(xref, file = "Res/20160803_kinaseSubstrateDatabaseParse/xref.uniprot.RDS")

if (!exists(".xref"))
  .xref <- readRDS("Res/20160803_kinaseSubstrateDatabaseParse/xref.uniprot.RDS")

accTo <- function(acc, to) {
  sub <- .xref[.xref$ACC %in% acc, ]
  sub <- sub[sub$source == to, ]
  sub <- unique(sub)
  id <- tapply(sub$value, sub$ACC, paste, collapse=";")
  id[acc]
}

gnToAcc <- function(gn) {
  sub <- .xref[.xref$value %in% gn, ]
  sub <- sub[sub$source == "Gene_Name", ]
  sub <- unique(sub)
  id <- tapply(sub$ACC, sub$value, paste, collapse=";")
  res <- id[gn]
  
  nagn <- gn[is.na(res)]
  sub2 <- .xref[.xref$value %in% nagn, ]
  sub2 <- sub2[sub2$source == "Gene_Synonym", ]
  id2 <- tapply(sub2$ACC, sub2$value, paste, collapse=";")
  res[is.na(res)] <- id2[nagn]
  
  res
}