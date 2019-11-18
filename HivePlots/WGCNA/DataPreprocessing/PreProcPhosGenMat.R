require(data.table)

annotation <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/crc65_psites_annotation.rds')
annotation_nci60_gluc <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/nci60_gluc_psites_annotation.rds')
annotation_fp_crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/crc65_fp_annotation.rds')
annotation_fp_nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/nci60_fp_annotation.rds')
annotation_fp_nci60_gluc <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/nci60_gluc_fp_annotation.rds')

dat <- readRDS("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Shiny/shinyImage0.3.RDS")
dat2 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20171218_prepareWgcnaMartin/dat.RDS')

allsites_crc65 <- dat2$crc65$mat

test <- as.matrix(annotation[rownames(allsites_crc65)[rownames(allsites_crc65)%in%annotation$rid],c(grep('_CRC65$', grep('Intensity ', colnames(annotation), value = T),value = T)),with=F])
test[test==0] <- NA
test <- log10(test)
rownames(test) <- annotation[rownames(allsites_crc65)[rownames(allsites_crc65)%in%annotation$rid],label]
colnames(test) <- gsub('Intensity ', '', colnames(test))
test <- test[,!colnames(test)%in%c('CoCM-1_CRC65')]
# diag(cor(test[rownames(test),],allsites_crc65[rownames(test),], use = 'p'))
test2 <- as.matrix(annotation_fp_crc65[rownames(allsites_crc65)[rownames(allsites_crc65)%in%annotation_fp_crc65$rid],c(grep('_CRC65$', grep('iBAQ ', colnames(annotation_fp_crc65), value = T),value = T)),with=F])
test2[test2==0] <- NA
test2 <- log10(test2)
rownames(test2) <- annotation_fp_crc65[rownames(allsites_crc65)[rownames(allsites_crc65)%in%annotation_fp_crc65$rid],firstgene]
colnames(test2) <- gsub('iBAQ ', '', colnames(test2))
test2 <- test2[,!colnames(test2)%in%c('CoCM-1_CRC65')]
# diag(cor(test2[rownames(test2),],allsites_crc65[rownames(test2),], use = 'p'))
allsites_crc65_nas <- rbind(test, test2[,colnames(test)])
# plot(allsites_crc65_nas[,"C10_CRC65"],allsites_crc65[,"C10_CRC65"])

# allsites_crc65_nas <- rbind(dat$tryp$fun_prepElnet(dat$tryp, panel = "CRC65", celllines = setdiff(colnames(dat$tryp$crc65_expr), "CoCM-1_CRC65"))$na_index,
#                             dat$fp.crc$fun_prepElnet(x = dat$fp.crc, celllines = setdiff(colnames(dat$tryp$crc65_expr), "CoCM-1_CRC65"))$na_index)
# allsites_crc65_nas <- allsites_crc65_nas[rownames(allsites_crc65),colnames(allsites_crc65)]
sanitycheck <- annotation_fp_crc65[id%in%setdiff(annotation_fp_crc65$id, as.integer(gsub(' .*', '', rownames(allsites_crc65)[1:which.min(diff(as.integer(gsub(' .*', '', rownames(allsites_crc65)))))]))),]
sanitycheck_small <- sanitycheck[,list(`Gene names`,`Only identified by site`,Reverse,`Potential contaminant`,`MS/MS count`,Intensity,Score, `Q-value`)]
rownames(allsites_crc65)[rownames(allsites_crc65)%in%annotation$rid] <- annotation[rownames(allsites_crc65)[rownames(allsites_crc65)%in%annotation$rid],label]
rownames(allsites_crc65)[rownames(allsites_crc65)%in%annotation_fp_crc65$rid] <- annotation_fp_crc65[rownames(allsites_crc65)[rownames(allsites_crc65)%in%annotation_fp_crc65$rid],firstgene]
allsites_crc65 <- allsites_crc65[intersect(rownames(allsites_crc65),c(annotation[`Gene names`!='',label],annotation_fp_crc65[`Gene names`!='',firstgene])),]
allsites_crc65_nas <- allsites_crc65_nas[rownames(allsites_crc65),colnames(allsites_crc65)]
diag(cor(allsites_crc65_nas,allsites_crc65, use = 'p'))
allsites_crc65_nas = !is.na(allsites_crc65_nas)
allsites_crc65_incnas <- allsites_crc65
saveRDS(allsites_crc65_incnas, file = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_crc65/RDS/allsites_nonas.rds', compress = T)
allsites_crc65_incnas[!allsites_crc65_nas] <- NA
saveRDS(allsites_crc65_incnas, file = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_crc65/RDS/allsites_incnas.rds', compress = T)
np_crc65 <- allsites_crc65_nas %*% t(allsites_crc65_nas)
saveRDS(np_crc65, file = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_crc65/RDS/nObs.rds', compress = T)
rm(np_crc65)

allsites_nci60_tryp <- dat2$nci60.tryp$mat

test <- as.matrix(annotation[rownames(allsites_nci60_tryp)[rownames(allsites_nci60_tryp)%in%annotation$rid],c(grep('_NCI60$', grep('Intensity ', colnames(annotation), value = T),value = T)),with=F])
test[test==0] <- NA
test <- log10(test)
rownames(test) <- annotation[rownames(allsites_nci60_tryp)[rownames(allsites_nci60_tryp)%in%annotation$rid],label]
colnames(test) <- gsub('Intensity ', '', colnames(test))
# test <- test[,!colnames(test)%in%c('CoCM-1_nci60_tryp')]
# diag(cor(test[rownames(test),],allsites_nci60_tryp[rownames(test),], use = 'p'))
test2 <- as.matrix(annotation_fp_nci60[rownames(allsites_nci60_tryp)[rownames(allsites_nci60_tryp)%in%annotation_fp_nci60$rid],c(grep('_NCI60$', grep('iBAQ ', colnames(annotation_fp_nci60), value = T),value = T)),with=F])
test2[test2==0] <- NA
test2 <- log10(test2)
rownames(test2) <- annotation_fp_nci60[rownames(allsites_nci60_tryp)[rownames(allsites_nci60_tryp)%in%annotation_fp_nci60$rid],firstgene]
colnames(test2) <- gsub('iBAQ ', '', colnames(test2))
# test2 <- test2[,!colnames(test2)%in%c('CoCM-1_nci60_tryp')]
# diag(cor(test2[rownames(test2),],allsites_nci60_tryp[rownames(test2),], use = 'p'))
allsites_nci60_tryp_nas <- rbind(test, test2[,colnames(test)])
# plot(allsites_nci60_tryp_nas[,"C10_nci60_tryp"],allsites_nci60_tryp[,"C10_nci60_tryp"])

# allsites_nci60_tryp_nas <- rbind(dat$tryp$fun_prepElnet(dat$tryp, panel = "NCI60", celllines = colnames(dat$tryp$nci60_expr))$na_index,
#                                  dat$fp.nci$fun_prepElnet(x = dat$fp.nci, celllines = colnames(dat$tryp$nci60_expr))$na_index)
# allsites_nci60_tryp_nas <- allsites_nci60_tryp_nas[rownames(allsites_nci60_tryp),colnames(allsites_nci60_tryp)]
rownames(allsites_nci60_tryp)[rownames(allsites_nci60_tryp)%in%annotation$rid] <- annotation[rownames(allsites_nci60_tryp)[rownames(allsites_nci60_tryp)%in%annotation$rid],label]
rownames(allsites_nci60_tryp)[rownames(allsites_nci60_tryp)%in%annotation_fp_nci60$rid] <- annotation_fp_nci60[rownames(allsites_nci60_tryp)[rownames(allsites_nci60_tryp)%in%annotation_fp_nci60$rid],firstgene]
allsites_nci60_tryp <- allsites_nci60_tryp[intersect(rownames(allsites_nci60_tryp),c(annotation[`Gene names`!='',label],annotation_fp_nci60[`Gene names`!='',firstgene])),]
allsites_nci60_tryp_nas <- allsites_nci60_tryp_nas[rownames(allsites_nci60_tryp),colnames(allsites_nci60_tryp)]
diag(cor(allsites_nci60_tryp_nas,allsites_nci60_tryp, use = 'p'))
allsites_nci60_tryp_nas = !is.na(allsites_nci60_tryp_nas)
allsites_nci60_tryp_incnas <- allsites_nci60_tryp
saveRDS(allsites_nci60_tryp_incnas, file = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_nci60_tryp/RDS/allsites_nonas.rds', compress = T)
allsites_nci60_tryp_incnas[!allsites_nci60_tryp_nas] <- NA
saveRDS(allsites_nci60_tryp_incnas, file = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_nci60_tryp/RDS/allsites_incnas.rds', compress = T)
np_nci60_tryp <- allsites_nci60_tryp_nas %*% t(allsites_nci60_tryp_nas)
saveRDS(np_nci60_tryp, file = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_nci60_tryp/RDS/nObs.rds', compress = T)
rm(np_nci60_tryp)

allsites_nci60_gluc <- dat2$nci60.gluc$mat

test <- as.matrix(annotation_nci60_gluc[rownames(allsites_nci60_gluc)[rownames(allsites_nci60_gluc)%in%annotation_nci60_gluc$rid],c(grep('___1$|___2$|___3$', grep('Intensity ', colnames(annotation_nci60_gluc), value = T),value = T, invert = T)),with=F])
test[test==0] <- NA
test <- log10(test)
rownames(test) <- annotation_nci60_gluc[rownames(allsites_nci60_gluc)[rownames(allsites_nci60_gluc)%in%annotation_nci60_gluc$rid],label]
colnames(test) <- paste0(gsub('Intensity ', '', colnames(test)), '_GluC')
test <- test[,!colnames(test)%in%c('Hela_GluC')]
colnames(test)[colnames(test)%in%c('LOX IMVI_GluC', 'NCIADRRES_GluC', 'UO31_GluC')] <- c('LOXIMVI_GluC', 'NCIADRES_GluC', 'U031_GluC')
# diag(cor(test[rownames(test),],allsites_nci60_gluc[rownames(test),], use = 'p'))
test2 <- as.matrix(annotation_fp_nci60_gluc[rownames(allsites_nci60_gluc)[rownames(allsites_nci60_gluc)%in%annotation_fp_nci60_gluc$rid],c(grep('___1$|___2$|___3$', grep('iBAQ ', colnames(annotation_fp_nci60_gluc), value = T),value = T, invert = T)),with=F])
test2[test2==0] <- NA
test2 <- log10(test2)
rownames(test2) <- annotation_fp_nci60_gluc[rownames(allsites_nci60_gluc)[rownames(allsites_nci60_gluc)%in%annotation_fp_nci60_gluc$rid],firstgene]
colnames(test2) <- gsub('iBAQ ', '', colnames(test2))
# test2 <- test2[,!colnames(test2)%in%c('CoCM-1_nci60_gluc')]
# diag(cor(test2[rownames(test2),],allsites_nci60_gluc[rownames(test2),], use = 'p'))
allsites_nci60_gluc_nas <- rbind(test, test2[,colnames(test)])
# plot(allsites_nci60_gluc_nas[,"C10_nci60_gluc"],allsites_nci60_gluc[,"C10_nci60_gluc"])

# allsites_nci60_gluc_nas <- rbind(dat$gluc$fun_prepElnet(dat$gluc, panel = "NCI60", celllines = colnames(dat$gluc$nci60_expr))$na_index,
#                                  dat$fp.nci$fun_prepElnet(x = dat$fp.nci, celllines = colnames(dat$gluc$nci60_expr))$na_index)
# allsites_nci60_gluc_nas <- allsites_nci60_gluc_nas[rownames(allsites_nci60_gluc),colnames(allsites_nci60_gluc)]
rownames(allsites_nci60_gluc)[rownames(allsites_nci60_gluc)%in%annotation_nci60_gluc$rid] <- annotation_nci60_gluc[rownames(allsites_nci60_gluc)[rownames(allsites_nci60_gluc)%in%annotation_nci60_gluc$rid],label]
rownames(allsites_nci60_gluc)[rownames(allsites_nci60_gluc)%in%annotation_fp_nci60_gluc$rid] <- annotation_fp_nci60_gluc[rownames(allsites_nci60_gluc)[rownames(allsites_nci60_gluc)%in%annotation_fp_nci60_gluc$rid],firstgene]
allsites_nci60_gluc <- allsites_nci60_gluc[intersect(rownames(allsites_nci60_gluc),c(annotation_nci60_gluc[`Gene names`!='',label],annotation_fp_nci60_gluc[`Gene names`!='',firstgene])),]
allsites_nci60_gluc_nas <- allsites_nci60_gluc_nas[rownames(allsites_nci60_gluc),colnames(allsites_nci60_gluc)]
diag(cor(allsites_nci60_gluc_nas,allsites_nci60_gluc, use = 'p'))
allsites_nci60_gluc_nas = !is.na(allsites_nci60_gluc_nas)
allsites_nci60_gluc_incnas <- allsites_nci60_gluc
saveRDS(allsites_nci60_gluc_incnas, file = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_nci60_gluc/RDS/allsites_nonas.rds', compress = T)
allsites_nci60_gluc_incnas[!allsites_nci60_gluc_nas] <- NA
saveRDS(allsites_nci60_gluc_incnas, file = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_nci60_gluc/RDS/allsites_incnas.rds', compress = T)
np_nci60_gluc <- allsites_nci60_gluc_nas %*% t(allsites_nci60_gluc_nas)
saveRDS(np_nci60_gluc, file = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_nci60_gluc/RDS/nObs.rds', compress = T)
rm(np_nci60_gluc)

# allsites_crc65 <- dat$tryp$fun_prepElnet(dat$tryp, panel = "CRC65", celllines = setdiff(colnames(dat$tryp$crc65_expr), "CoCM-1_CRC65"))
# regsites_crc65 <- dat$tryp$fun_prepElnet(dat$tryp, panel = "CRC65", celllines = setdiff(colnames(dat$tryp$crc65_expr), "CoCM-1_CRC65"), regsitesOnly = T)
# fp_crc65 <- dat$fp.crc$fun_prepElnet(dat$fp.crc, celllines = setdiff(colnames(dat$fp.crc$ibaq), "CoCM-1_CRC65"))
# 
# allsites_crc65_NAs <- allsites_crc65$imputed
# rownames(allsites_crc65_NAs) <- annotation[rownames(allsites_crc65_NAs),label]
# allsites_crc65_NAs <- allsites_crc65_NAs[intersect(rownames(allsites_crc65_NAs),
#                                                    annotation[`Gene names`!='',label]),]
# 
# regsites_crc65_NAs <- regsites_crc65$imputed
# # regsites_NAs[regsites$na_index] <- NA
# rownames(regsites_crc65_NAs) <- annotation[rownames(regsites_crc65_NAs),label]
# regsites_crc65_NAs <- regsites_crc65_NAs[intersect(rownames(regsites_crc65_NAs),
#                                                    annotation[`Gene names`!='',label]),]
# 
# fp_crc65_NAs <- fp_crc65$imputed
# rownames(fp_crc65_NAs) <- annotation_fp_crc65[rownames(fp_crc65_NAs),firstgene]
# fp_crc65_NAs <- fp_crc65_NAs[intersect(rownames(fp_crc65_NAs),
#                                        annotation_fp_crc65[`Gene names`!='',firstgene]),]
# 
# 
# allsites_nci60 <- dat$tryp$fun_prepElnet(dat$tryp, panel = "NCI60", celllines = colnames(dat$tryp$nci60_expr))
# regsites_nci60 <- dat$tryp$fun_prepElnet(dat$tryp, panel = "NCI60", celllines = colnames(dat$tryp$nci60_expr), regsitesOnly = T)
# fp_nci60 <- dat$fp.nci$fun_prepElnet(dat$fp.nci, celllines = colnames(dat$fp.nci$ibaq))
# 
# allsites_nci60_NAs <- allsites_nci60$imputed
# rownames(allsites_nci60_NAs) <- annotation[rownames(allsites_nci60_NAs),label]
# allsites_nci60_NAs <- allsites_nci60_NAs[intersect(rownames(allsites_nci60_NAs),
#                                                    annotation[`Gene names`!='',label]),]
# 
# regsites_nci60_NAs <- regsites_nci60$imputed
# # regsites_NAs[regsites$na_index] <- NA
# rownames(regsites_nci60_NAs) <- annotation[rownames(regsites_nci60_NAs),label]
# regsites_nci60_NAs <- regsites_nci60_NAs[intersect(rownames(regsites_nci60_NAs),
#                                                    annotation[`Gene names`!='',label]),]
# 
# fp_nci60_NAs <- fp_nci60$imputed
# rownames(fp_nci60_NAs) <- annotation_fp_nci60[rownames(fp_nci60_NAs),firstgene]
# fp_nci60_NAs <- fp_nci60_NAs[intersect(rownames(fp_nci60_NAs),
#                                        annotation_fp_nci60[`Gene names`!='',firstgene]),]
# 
# 
# allsites_nci60_gluc <- dat$gluc$fun_prepElnet(dat$gluc, celllines = colnames(dat$gluc$nci60_expr))
# regsites_nci60_gluc <- dat$gluc$fun_prepElnet(dat$gluc, celllines = colnames(dat$gluc$nci60_expr), regsitesOnly = T)
# fp_nci60_gluc <- dat$fp.nci$fun_prepElnet(dat$fp.nci, celllines = colnames(dat$fp.nci$ibaq))
# 
# allsites_nci60_gluc_NAs <- allsites_nci60_gluc$imputed
# rownames(allsites_nci60_gluc_NAs) <- annotation[rownames(allsites_nci60_gluc_NAs),label]
# allsites_nci60_gluc_NAs <- allsites_nci60_gluc_NAs[intersect(rownames(allsites_nci60_gluc_NAs),
#                                                    annotation[`Gene names`!='',label]),]
# 
# regsites_nci60_gluc_NAs <- regsites_nci60_gluc$imputed
# # regsites_NAs[regsites$na_index] <- NA
# rownames(regsites_nci60_gluc_NAs) <- annotation[rownames(regsites_nci60_gluc_NAs),label]
# regsites_nci60_gluc_NAs <- regsites_nci60_gluc_NAs[intersect(rownames(regsites_nci60_gluc_NAs),
#                                                    annotation[`Gene names`!='',label]),]
# 
# fp_nci60_gluc_NAs <- fp_nci60_gluc$imputed
# rownames(fp_nci60_gluc_NAs) <- annotation_fp_nci60_gluc[rownames(fp_nci60_gluc_NAs),firstgene]
# fp_nci60_gluc_NAs <- fp_nci60_gluc_NAs[intersect(rownames(fp_nci60_gluc_NAs),
#                                        annotation_fp_nci60_gluc[`Gene names`!='',firstgene]),]


# pkinfam1 <- fread("/media/msdata5/users_files/martin/phd/data/ms/24_CRC64_1and2and3_120ppm_1.4.1.2_peakproperties/analysis/pkinfam1.txt")
# setkey(pkinfam1, Name)
# pkinfam1 <- pkinfam1[!""]
# 
# fp_NAs_kin <- fp_NAs[rownames(fp_NAs)%in%pkinfam1$Name,]

traits <- dat$tryp$crc65_xref
mutnames <- names(which(colSums(!traits[,seq(match('De.Sousa.E.Melo.subtype', colnames(traits))+1,ncol(traits),1)], na.rm = T)!=0))

traits1 <- traits[,c("CMS.mRNA",
                     "Marisa.subtype",
                     "Budinska.subtype",
                     "Roepman.subtype",
                     "Sadanandam.subtype",
                     "De.Sousa.E.Melo.subtype")]
# for (i in colnames(traits1)) {
#   traits1[,i] <- factor(traits1[,i],
#                         labels = unique(traits1[,i])[!is.na(unique(traits1[,i]))],
#                         levels = unique(traits1[,i])[!is.na(unique(traits1[,i]))])
# }
# traits1 <- apply(traits1, 2, function(i) factor(i, labels = unique(i)[!is.na(unique(i))], levels = unique(i)[!is.na(unique(i))]))
traits1 <- Reduce(cbind, lapply(colnames(traits1), function(i) {
  nm <- gsub('\\.mRNA|\\.subtype', '', i)
  res <- as.data.frame.matrix(table(rownames(traits1), traits1[,i]))
  res[rowSums(res)==0,] <- NA
  colnames(res) <- paste0(nm, '_', colnames(res))
  return(res)
}
)
)

# rownames(traits1) <- rownames(traits)

traits2 <- traits[,mutnames]
traits2 <- apply(traits2, 2, function(i) as.integer(i))
rownames(traits2) <- rownames(traits)

traits <- cbind(traits1, traits2)

# traits <- traits[rownames(traits)%in%colnames(allsites_crc65),]
colnames(traits) <- gsub('\\.', ' ', colnames(traits))
saveRDS(traits, '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/traits_crc65.rds', compress = T)


traits_nci60 <- dat2$nci60.tryp$xref
# invisible(lapply(colnames(traits_nci60), function(i) traits_nci60[eval(as.name(i))=='',eval(i):=NA]))
traits_nci60[traits_nci60==''] <- NA
traits_nci60 <- data.frame(traits_nci60, check.names = F)
factors <- names(which(sapply(colnames(traits_nci60), function(i) class(traits_nci60[,i]))=='character'))
traits1_nci60 <- traits_nci60[,factors]

traits1_nci60 <- Reduce(cbind, lapply(colnames(traits1_nci60), function(i) {
  nm <- gsub('\\.mRNA|\\.subtype', '', i)
  res <- as.data.frame.matrix(table(rownames(traits1_nci60), traits1_nci60[,i]))
  res[rowSums(res)==0,] <- NA
  colnames(res) <- paste0(nm, '_', colnames(res))
  return(res)
}
)
)

traits_nci60 <- cbind(traits1_nci60, traits_nci60[,!colnames(traits_nci60)%in%factors])
saveRDS(traits_nci60, '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/traits_nci60.rds', compress = T)

rm(dat)
rm(dat2)
gc()