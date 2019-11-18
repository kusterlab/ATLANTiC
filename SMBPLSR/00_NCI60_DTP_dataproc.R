library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)
library(omic3plus)
library(impute)

# source("../R/elnetFit.R")
source("R/drugResponseProcess/extract.drugmat.nci60.R")


dim(drug.nci60$dtp.gi50)

prot.gluc <- readRDS("Res/20170719_nci60.proteinGroups.gluc/nci60.proteinGroups.gluc.RDS")
prot.tryp <- readRDS("Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
site.gluc <- readRDS("Res/20170110_mq15processLOD2/sites.gluC.RDS")
site.tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
#

ii <- setdiff(intersect(colnames(drug.nci60$dtp.gi50), colnames(prot.tryp$intensity)), "MDAMB468_NCI60")

pg <- prot.gluc$fun_prepElnet(prot.gluc, celllines = gsub("NCI60", "GluC", ii))
pt <- prot.tryp$fun_prepElnet(prot.tryp, celllines = ii)
sg <- site.gluc$fun_prepElnet(site.gluc, celllines = gsub("_NCI60", "", ii))
st <- site.tryp$fun_prepElnet(site.tryp, celllines = ii, panel = "NCI60")

dtp <- drug.nci60$dtp.gi50[, ii]

# impute missing value
dtpi <- impute.knn(dtp)
dim(dtpi$data)

dat <- list( y = dtpi, x = list(prot_gluc = pg, prot_tryp = pt, 
                               site_gluc = sg, site.tryp = st))

dat$drugannot = drug.nci60$dtp.annot
dat$xref = site.tryp$nci60_xref[colnames(dat$x$site.tryp$imputed), ]

saveRDS(dat, file = "Res/20170904_concordanceData/nci60_dtp.RDS")



