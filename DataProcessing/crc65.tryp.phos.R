tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")

tryp.id <- paste(sapply(strsplit(tryp$annot$`Leading proteins`, "-|;"), "[", 1), 
                 tryp$annot$Position, 
                 sep = "_")

xref <- tryp$crc65_xref
tryp.list <- lapply(xref$name_MQ, function(x) {
  tryp.id[!is.na(tryp$crc65_expr[, x])]
})

res <- list(cellline = structure(sapply(tryp.list, length), names = xref$name_MQ), 
            tissue = rbind(nsample = nrow(xref),
                           tisid = length(unique(unlist(tryp.id)))))

saveRDS(res, file = "Res/20170721_numbersFinalPlotData/numberid_crc65_PHOS.RDS")


