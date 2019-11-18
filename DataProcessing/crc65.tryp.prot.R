# ======================================== FP ================================================
### crc65 FP
tryp <- readRDS(file = "Res/20170110_mq15processLOD2/crc65.proteinGroups.RDS")
tryp.id <- sapply(strsplit(tryp$annot$`Protein IDs`, ";|-"), "[", 1)

xref <- tryp$xref
tryp.list <- lapply(xref$name_MQ, function(x) {
  tryp.id[!is.na(tryp$intensity[, x])]
})

crc65.fp <- tryp.list
length(unique(unlist(crc65.fp)))


###  NCI60 FP
tryp <- readRDS(file = "Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
gluc <- readRDS(file = "Res/20170719_nci60.proteinGroups.gluc/nci60.proteinGroups.gluc.RDS")

tryp.id <- sapply(strsplit(tryp$annot$`Protein IDs`, ";|-"), "[", 1)
gluc.id <- sapply(strsplit(gluc$annot$`Protein IDs`, ";|-"), "[", 1)

#
xref <- tryp$xref
xref <- xref[order(xref$too.short), ]

#### count
tryp.list <- lapply(xref$name_MQ, function(x) {
  tryp.id[!is.na(tryp$intensity[, x])]
})
gluc.list <- lapply(gsub("_NCI60", "_GluC", xref$name_MQ), function(x) {
  gluc.id[!is.na(gluc$intensity[, x])]
})

nci.tryp.fp <- tryp.list
nci.gluc.fp <- gluc.list

## total number of proteins
length(unique(c(unlist(nci.tryp.fp), unlist(nci.gluc.fp), unlist(crc65.fp))))


