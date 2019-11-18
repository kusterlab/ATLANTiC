library(reshape2)
###
corCRC65 <- readRDS("/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/res/20180302_wgcna_allsites_allproteins_incnas_crc65_nonas_crc65/RDS/corandp.rds")
corCRC65 <- corCRC65$NoImputation
gc()

cormat <- corCRC65$cor
cormat[corCRC65$nObs<7] <- NA
rm(corCRC65)
gc()
cormat[ abs(cormat) < 0.5 ] <- NA
gc()
cormat[ lower.tri(cormat, diag = TRUE) ] <- NA
gc()
gc()
cormat <- melt(cormat, na.rm = TRUE)
gc()

saveRDS(cormat, file = "Res/20180314_wgcnaSupTables/crc65_cormat.RDS")
rm(cormat)
gc()

###
corNCI60 <- readRDS("/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/res/20180308_wgcna_allsites_allproteins_incnas_nci60_nonas_nci60/RDS/corandp.rds")
corNCI60 <- corNCI60$NoImputation
gc()

cormat <- corNCI60$cor
gc()
cormat[corNCI60$nObs<7] <- NA
rm(corNCI60)
gc()
cormat[ abs(cormat) < 0.5 ] <- NA
gc()
cormat[ lower.tri(cormat, diag = TRUE) ] <- NA
gc()
gc()
gc()
cormat <- melt(cormat, na.rm = TRUE)
gc()
saveRDS(cormat, file = "Res/20180314_wgcnaSupTables/nci60_cormat.RDS")