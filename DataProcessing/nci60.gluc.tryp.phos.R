tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
gluc <- readRDS("Res/20170110_mq15processLOD2/sites.gluC.RDS")


#
xref <- gluc$nci60_xref
xref <- xref[order(xref$too.short), ]

#
gluc.id <- paste(sapply(strsplit(gluc$annot$`Leading proteins`, "-|;"), "[", 1), 
                 gluc$annot$Position, 
                 sep = "_")
tryp.id <- paste(sapply(strsplit(tryp$annot$`Leading proteins`, "-|;"), "[", 1), 
                 tryp$annot$Position, 
                 sep = "_")

#### list
tryp.list <- lapply(xref$name_MQ, function(x) {
  tryp.id[!is.na(tryp$nci60_expr[, x])]
})
gluc.list <- lapply(gsub("_NCI60", "", xref$name_MQ), function(x) {
  gluc.id[!is.na(gluc$nci60_expr[, x])]
})


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
saveRDS(res, file = "Res/20170721_numbersFinalPlotData/numberid_nci60_PHOS.RDS")


x <- readRDS("Res/20170721_numbersFinalPlotData/numberid_nci60_PHOS.RDS")
x <- readRDS("Res/20170721_numbersFinalPlotData/numberid_nci60_FP.RDS")

x
