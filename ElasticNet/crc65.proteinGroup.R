library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)

source("../R/elnetFit.R")
source("R/drugResponseProcess/extract.drugmat.crc65.R")
slist <- readRDS("Res/20170110_mq15processLOD2/crc65.proteinGroups.RDS")

colnames(drug.crc65$CTRP)
colnames(slist$ibaq)

runMod <- function( drugMat, sites, regsitesOnly = FALSE) {# the drug matrix, rows are drugs columns are samples
  
  if (! identical(colnames(drugMat), gsub("_CRC65", "", colnames(sites$ibaq))))
    stop("check columns order.")
  
  fit <- list()
  for (i in 1:nrow(drugMat)) {
    print(i)
    gi50 <- drugMat[i, ]
    iSample <- !is.na(gi50)
    
    if (sum(iSample) < 12) {
      fit[[i]] <- NA
      next()
    }
    
    ixmat <- sites$fun_prepElnet(x = sites, celllines = iSample)
    fit[[i]] <- try(elnetFit(t(ixmat$imputed), 
                             y = gi50[iSample], 
                             standardizex = TRUE,
                             nboot = 100, 
                             permpvalue = 0, 
                             x.na.index = t(ixmat$na_index),
                             times=1, 
                             cores = 22, 
                             write.table = FALSE), silent = TRUE)
  }
  names(fit) <- rownames(drugMat)
  list(fit=fit, mat=ixmat)
}

## all phospho-sites
ccle <- runMod(drugMat = drug.crc65$CCLE, sites = slist)
saveRDS(ccle, file = "Res/20170110_mq15processLOD2/crc65.ccle.proteinGroups.RDS")

medico <- runMod(drugMat = drug.crc65$MEDICO, sites = slist)
saveRDS(medico, file = "Res/20170110_mq15processLOD2/crc65.medico.proteinGroups.RDS")

gdsc <- runMod(drugMat = drug.crc65$GDSC, sites = slist)
saveRDS(gdsc, file = "Res/20170110_mq15processLOD2/crc65.gdsc.proteinGroups.RDS")

ctrp <- runMod(drugMat = drug.crc65$CTRP, sites = slist)
saveRDS(ctrp, file = "Res/20170110_mq15processLOD2/crc65.ctrp.proteinGroups.RDS")

