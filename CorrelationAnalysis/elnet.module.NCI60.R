library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)
library(psych)
library(reshape2)

source("../R/elnetFit.R")
source("R/drugResponseProcess/extract.drugmat.nci60.R")
wgcna.pc1 <- readRDS("Res/20170727_wgcna/summary.pc1.RDS")

m.nci60 <- c()
mat <- t(sapply(wgcna.pc1$site.nci60.tryp, "[[", "pc1"))
m.nci60 <- rbind(m.nci60, mat)
mat <- t(sapply(wgcna.pc1$site.nci60.gluc, "[[", "pc1"))
m.nci60 <- rbind(m.nci60, mat)
mat <- t(sapply(wgcna.pc1$fp.nci60.tryp, "[[", "pc1"))
m.nci60 <- rbind(m.nci60, mat)
mat <- t(sapply(wgcna.pc1$fp.nci60.gluc, "[[", "pc1"))
m.nci60 <- rbind(m.nci60, mat)
dim(m.nci60)

identical(colnames(m.nci60), colnames(drug.nci60$dtp.gi50))

## Define function

runMod <- function( drugMat, sites, regsitesOnly = FALSE) {# the drug matrix, rows are drugs columns are samples
  
  if (! identical(colnames(drugMat), colnames(sites)))
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
    
    xmat <- sites[, iSample]
    fit[[i]] <- try(elnetFit(t(xmat), 
                             y = gi50[iSample], 
                             standardizex = TRUE,
                             nboot = 100, 
                             permpvalue = 0, 
                             x.na.index = matrix(FALSE, nrow(xmat), ncol(xmat)),
                             times=1, 
                             cores = 22, 
                             write.table = FALSE), silent = TRUE)
  }
  names(fit) <- rownames(drugMat)
  list(fit=fit, mat=xmat)
}


## all phospho-sites
ccle <- runMod(drugMat = drug.nci60$ccle.auc, sites = m.nci60)
saveRDS(ccle, file = "Res/20170731_weight.comboScore/elnet.allModules.nci60_ccle.RDS")

dtp <- runMod(drugMat = drug.nci60$dtp.gi50, sites = m.nci60)
saveRDS(dtp, file = "Res/20170731_weight.comboScore/elnet.allModules.nci60_dtp.RDS")

gdsc <- runMod(drugMat = drug.nci60$gdsc.auc, sites = m.nci60)
saveRDS(gdsc, file = "Res/20170731_weight.comboScore/elnet.allModules.nci60_gdsc.RDS")

ctrp <- runMod(drugMat = drug.nci60$ctrp.auc, sites = m.nci60)
saveRDS(ctrp, file = "Res/20170731_weight.comboScore/elnet.allModules.nci60_ctrp.RDS")
