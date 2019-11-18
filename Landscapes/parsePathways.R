require(data.table)
require(reshape2)
require(scales)

geoMean <- function (x, na.rm = TRUE) 
{
  if (is.null(nrow(x))) {
    exp(mean(log(x), na.rm = TRUE))
  }
  else {
    exp(apply(log(x), 2, mean, na.rm = na.rm))
  }
}
adjustScores <- function(scorelist, scores = c('Kinases (by Kinase abundance)',
                                               'Kinases (by Kinase phosphorylation)',
                                               'Kinases (by Substrate phosphorylation)'), adjust = F, mediancenter = F) {
  overlap <- Reduce(intersect, lapply(scorelist[scores], rownames))
  cmbl <- lapply(scorelist[scores], function(i) i[overlap,])
  if (mediancenter) {
    cmbl <- lapply(cmbl, mediancenter)
  }
  if (adjust) {
    require(sva)
    mgd <- do.call(cbind, cmbl)
    btc <- unlist(lapply(seq_along(scores), function(i) rep(i, ncol(cmbl[[scores[i]]]))))
    cmb <- ComBat(mgd, btc)
    cmbl <- lapply(unique(btc), function(i) cmb[,btc==i])
    names(cmbl) <- scores
  }
  return(cmbl)
}
mediancenter <- function(x, refrows=NULL) {
  if (!is.null(refrows)) {
    gm <- apply(x[refrows,], 2, median, na.rm = T)
    x <- sweep(x, 2, gm, `-`)+median(x[refrows,], na.rm=T)
  } else {
    gm <- apply(x, 2, median, na.rm=T)
    x <- sweep(x, 2, gm, `-`)+median(x, na.rm=T)
  }
  return(x)
}
pkinfam1 <- readRDS('/media/kusterlab/internal_projects/active/CRC65/data/24_CRC64_1and2and3_120ppm_1.4.1.2_peakproperties/analysis/pkinfam1.rda')
setkey(pkinfam1, Name)
pkinfam1 <- pkinfam1[!""]

nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist_nci60.RDS') # from zscore.R
# nci602 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist_nci60_sub.RDS')
# nci60$integ <- nci602$integ
names(nci60) <- c('NetworKIN', 'Pathways', 'Kinases (by Substrate phosphorylation)', 'SIGNOR', 'Integrated_NetworKIN')

ev_phosprots_nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/ev_phosprots_nci60.rds')
phoskin_nci60 <- ev_phosprots_nci60[rownames(ev_phosprots_nci60)%in%pkinfam1$Name,]
for (i in 1:nrow(phoskin_nci60)) {
  phoskin_nci60[i,][is.na(phoskin_nci60[i,])] <- min(phoskin_nci60[i,], na.rm = T) - log10(2)
}
phoskin_nci60 <- t(scale(t(phoskin_nci60)))
nci60[['Kinases (by Kinase phosphorylation)']] <- phoskin_nci60

ev_prots_nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/ev_prots_nci60.rds')
protkin_nci60 <- ev_prots_nci60[rownames(ev_prots_nci60)%in%pkinfam1$Name,]
for (i in 1:nrow(protkin_nci60)) {
  protkin_nci60[i,][is.na(protkin_nci60[i,])] <- min(protkin_nci60[i,], na.rm = T) - log10(2)
}
protkin_nci60 <- t(scale(t(protkin_nci60)))
colnames(protkin_nci60) <- gsub('_NCI60', '', colnames(protkin_nci60))
nci60[['Kinases (by Kinase abundance)']] <- protkin_nci60

for (i in seq_along(nci60)) {
  colnames(nci60[[i]]) <- gsub('_NCI60', '', colnames(nci60[[i]]))
}

comb <- adjustScores(nci60, adjust = F, mediancenter = F)

# overlap <- Reduce(intersect, list(rownames(nci60$`Kinases (by Kinase abundance)`),
#                                   rownames(nci60$`Kinases (by Kinase phosphorylation)`),
#                                   rownames(nci60$`Kinases (by Substrate phosphorylation)`)))

abu <- setDT(melt(comb$`Kinases (by Kinase abundance)`))
setnames(abu, c('Kinase', 'Cellline', 'Protein'))
pho <- setDT(melt(comb$`Kinases (by Kinase phosphorylation)`))
setnames(pho, c('Kinase', 'Cellline', 'Phosphoprotein'))
sub <- setDT(melt(comb$`Kinases (by Substrate phosphorylation)`))
setnames(sub, c('Kinase', 'Cellline', 'Substrates'))
merged <- merge(merge(abu, pho, by = c('Kinase', 'Cellline')), sub, by = c('Kinase', 'Cellline'))
# ATLANTIC <- acast(merged[,list(ATLANTIC=mean(c(Protein,Phosphoprotein,Substrates))),by=list(Kinase, Cellline)], Kinase~Cellline, value.var = 'ATLANTIC')
merged[,Protein:=rescale(Protein, c(0,1))]
merged[,Phosphoprotein:=rescale(Phosphoprotein, c(0,1))]
merged[,Substrates:=rescale(Substrates, c(0,1))]
ATLANTIC <- acast(merged[,list(ATLANTIC=(Protein*Phosphoprotein*Substrates)^(1/3)),by=list(Kinase, Cellline)], Kinase~Cellline, value.var = 'ATLANTIC')
nci60[['Kinases (by ATLANTIC Score)']] <- ATLANTIC
nci60[c('NetworKIN', 'SIGNOR', 'Integrated_NetworKIN')] <- NULL

crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist_crc65.RDS') # from zscore.R
# crc652 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist_crc65_sub.RDS')
# crc65$integ <- crc652$integ

names(crc65) <- c('NetworKIN', 'Pathways', 'Kinases (by Substrate phosphorylation)', 'SIGNOR', 'Integrated_NetworKIN')

ev_phosprots_crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/ev_phosprots_crc65.rds')
phoskin_crc65 <- ev_phosprots_crc65[rownames(ev_phosprots_crc65)%in%pkinfam1$Name,]
for (i in 1:nrow(phoskin_crc65)) {
  phoskin_crc65[i,][is.na(phoskin_crc65[i,])] <- min(phoskin_crc65[i,], na.rm = T) - log10(2)
}
phoskin_crc65 <- t(scale(t(phoskin_crc65)))
colnames(phoskin_crc65) <- gsub('_CRC65', '', colnames(phoskin_crc65))
crc65[['Kinases (by Kinase phosphorylation)']] <- phoskin_crc65

ev_prots_crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/ev_prots_crc65.rds')
protkin_crc65 <- ev_prots_crc65[rownames(ev_prots_crc65)%in%pkinfam1$Name,]
for (i in 1:nrow(protkin_crc65)) {
  protkin_crc65[i,][is.na(protkin_crc65[i,])] <- min(protkin_crc65[i,], na.rm = T) - log10(2)
}
protkin_crc65 <- t(scale(t(protkin_crc65)))
crc65[['Kinases (by Kinase abundance)']] <- protkin_crc65

for (i in seq_along(crc65)) {
  colnames(crc65[[i]]) <- gsub('_CRC65', '', colnames(crc65[[i]]))
}

comb <- adjustScores(crc65, adjust = F, mediancenter = F)

# overlap <- Reduce(intersect, list(rownames(crc65$`Kinases (by Kinase abundance)`),
#                                   rownames(crc65$`Kinases (by Kinase phosphorylation)`),
#                                   rownames(crc65$`Kinases (by Substrate phosphorylation)`)))

abu <- setDT(melt(comb$`Kinases (by Kinase abundance)`))
setnames(abu, c('Kinase', 'Cellline', 'Protein'))
pho <- setDT(melt(comb$`Kinases (by Kinase phosphorylation)`))
setnames(pho, c('Kinase', 'Cellline', 'Phosphoprotein'))
sub <- setDT(melt(comb$`Kinases (by Substrate phosphorylation)`))
setnames(sub, c('Kinase', 'Cellline', 'Substrates'))
merged <- merge(merge(abu, pho, by = c('Kinase', 'Cellline')), sub, by = c('Kinase', 'Cellline'))
# ATLANTIC <- acast(merged[,list(ATLANTIC=mean(c(Protein, Phosphoprotein, Substrates))),by=list(Kinase, Cellline)], Kinase~Cellline, value.var = 'ATLANTIC')
merged[,Protein:=rescale(Protein, c(0,1))]
merged[,Phosphoprotein:=rescale(Phosphoprotein, c(0,1))]
merged[,Substrates:=rescale(Substrates, c(0,1))]
ATLANTIC <- acast(merged[,list(ATLANTIC=(Protein*Phosphoprotein*Substrates)^(1/3)),by=list(Kinase, Cellline)], Kinase~Cellline, value.var = 'ATLANTIC')
crc65[['Kinases (by ATLANTIC Score)']] <- ATLANTIC
crc65[c('NetworKIN', 'SIGNOR', 'Integrated_NetworKIN')] <- NULL

path <- list(CRC65=crc65, NCI60=nci60)
names(path$CRC65)[names(path$CRC65)=="Kinases (by ATLANTIC Score)"] <- "Kinases (by ATLANTiC Score)"
names(path$NCI60)[names(path$NCI60)=="Kinases (by ATLANTIC Score)"] <- "Kinases (by ATLANTiC Score)"
for (i in names(path)) {
  for (j in names(path[[i]])) {
    path[[i]][[j]] <- rescale(path[[i]][[j]], c(0,1))
  }
}
saveRDS(path, '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/zscore_data_new.RDS', compress = T)

