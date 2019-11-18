library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)

source("../R/elnetFit.R")
source("R/drugResponseProcess/extract.drugmat.crc65.R")
slist <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")


###

runMod <- function( drugMat, sites, regsitesOnly = FALSE) {
  
  if (! identical(colnames(drugMat), gsub("_CRC65", "", colnames(sites$crc65_expr))) )
    stop("check columns order.")
  
  fit <- list()
  for (i in 1:nrow(drugMat)) {
    print(i)
    gi50 <- drugMat[i, ]
    iSample <- !is.na(gi50)
    iSample[colnames(drugMat) == "CoCM-1"] <- FALSE # remove abnormal sample
    
    if (sum(iSample) < 12) {
      fit[[i]] <- NA
      next()
    }
    
    ixmat <- sites$fun_prepElnet(x = sites, panel = "CRC65", celllines = iSample, 
                                 regsitesOnly = regsitesOnly)
    
    fit[[i]] <- try( elnetFit(t(ixmat$imputed), 
                              y = gi50[iSample], 
                              nboot = 100,
                              standardizex = TRUE,
                              permpvalue = 0, 
                              x.na.index = t(ixmat$na_index),
                              times=1, 
                              cores = 20, 
                              write.table = FALSE) )
  }
  names(fit) <- rownames(drugMat)
  list(fit=fit, mat=NA)
}

## all phospho-sites
ccle <- runMod(drugMat = drug.crc65$CCLE, sites = slist)
saveRDS(ccle, file = "Res/20170110_mq15processLOD2/crc65.ccle.allsites.RDS")

gdsc <- runMod(drugMat = drug.crc65$GDSC, sites = slist)
saveRDS(gdsc, file = "Res/20170110_mq15processLOD2/crc65.gdsc.allsites.RDS")

ctrp <- runMod(drugMat = drug.crc65$CTRP, sites = slist)
saveRDS(ctrp, file = "Res/20170110_mq15processLOD2/crc65.ctrp.allsites.RDS")

medico <- runMod(drugMat = drug.crc65$MEDICO, sites = slist)
saveRDS(medico, file = "Res/20170110_mq15processLOD2/crc65.medico.allsites.RDS")


## regulatory sites only
ccle.rs <- runMod(drugMat = drug.crc65$CCLE, sites = slist, regsitesOnly = TRUE) 
saveRDS(ccle.rs, file = "Res/20170110_mq15processLOD2/crc65.ccle.regsites.RDS")

gdsc.rs <- runMod(drugMat = drug.crc65$GDSC, sites = slist, regsitesOnly = TRUE)
saveRDS(gdsc.rs, file = "Res/20170110_mq15processLOD2/crc65.gdsc.regsites.RDS")

ctrp.rs <- runMod(drugMat = drug.crc65$CTRP, sites = slist, regsitesOnly = TRUE)
saveRDS(ctrp.rs, file = "Res/20170110_mq15processLOD2/crc65.ctrp.regsites.RDS")

medico.rs <- runMod(drugMat = drug.crc65$MEDICO, sites = slist, regsitesOnly = TRUE)
saveRDS(medico.rs, file = "Res/20170110_mq15processLOD2/crc65.medico.regsites.RDS")


