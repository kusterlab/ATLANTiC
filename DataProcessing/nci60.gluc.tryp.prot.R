tryp <- readRDS(file = "Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
gluc <- readRDS(file = "Res/20170719_nci60.proteinGroups.gluc/nci60.proteinGroups.gluc.RDS")

tryp.id <- sapply(strsplit(tryp$annot$`Protein IDs`, ";|-"), "[", 1)
gluc.id <- sapply(strsplit(gluc$annot$`Protein IDs`, ";|-"), "[", 1)

#
xref <- tryp$xref
xref <- xref[order(xref$too.short), ]

#### list
tryp.list <- lapply(xref$name_MQ, function(x) {
  tryp.id[!is.na(tryp$intensity[, x])]
})
gluc.list <- lapply(gsub("_NCI60", "_GluC", xref$name_MQ), function(x) {
  gluc.id[!is.na(gluc$intensity[, x])]
})


length(unique(c(unlist(tryp.list), unlist(gluc.list))))

## common ids
both <- mapply(function(x, y) {
  length(intersect(x, y))
}, x = tryp.list, y = gluc.list)

trypOnly <- mapply(function(x, y) {
  length(setdiff(x, y))
}, x = tryp.list, y = gluc.list)

glucOnly <- mapply(function(x, y) {
  length(setdiff(y, x))
}, x = tryp.list, y = gluc.list)


# tissue specific
too <- unique(xref$too.short)
tisid <- sapply(too, function(x) {
  i <- xref$too.short == x
  length(unique(c(unlist(gluc.list[i]), unlist(tryp.list[i]))))
})

nsample <- table(xref$too.short)[too]

tis <- rbind(nsample, tisid)
cl <- rbind(both, trypOnly, glucOnly)
colnames(cl) <- xref$name_MQ

res <- list(cellline = cl, tissue = tis)
saveRDS(res, file = "Res/20170721_numbersFinalPlotData/numberid_nci60_FP.RDS")


