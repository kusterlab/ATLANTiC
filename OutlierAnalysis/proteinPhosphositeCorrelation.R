slist <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
prot <- readRDS("Res/20170110_mq15processLOD2/crc65.proteinGroups.RDS")

dd <- slist$crc65_expr[which(slist$annot$`Gene names` == "ALK"), ]
vv <- prot$ibaq[prot$annot$`Gene name` == "ALK", ]
par(mar = c(10, 4, 1, 1))
barplot(sort(10^vv, na.last = TRUE), las = 2)



dd <- slist$crc65_expr[which(slist$annot$`Gene names` == "EML4"), ]
vv <- prot$ibaq[prot$annot$`Gene name` == "EML4", ]
par(mar = c(10, 4, 1, 1))
barplot(sort(colSums(10^dd, na.rm = TRUE), na.last = TRUE), las = 2)
barplot(sort(10^vv, na.last = TRUE), las = 2)
