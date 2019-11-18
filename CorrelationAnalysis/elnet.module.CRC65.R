library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)
library(psych)
library(reshape2)

source("../R/elnetFit.R")
source("R/drugResponseProcess/extract.drugmat.crc65.R")
wgcna.pc1 <- readRDS("Res/20170727_wgcna/summary.pc1.RDS")


# CRC 65
m.crc65 <- c()
mat <- t(sapply(wgcna.pc1$site.crc65.tryp, "[[", "pc1"))
m.crc65 <- rbind(m.crc65, mat)
mat <- t(sapply(wgcna.pc1$fp.crc65.tryp, "[[", "pc1"))
m.crc65 <- rbind(m.crc65, mat)
dim(m.crc65)
colnames(m.crc65) <- gsub("_CRC65", "", colnames(m.crc65))

runMod <- function( drugMat, sites, regsitesOnly = FALSE) {
  
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

ii <- intersect(colnames(drug.crc65$CCLE), colnames(m.crc65))


## all phospho-sites
ll <- list()
ll$ccle <- runMod(drugMat = drug.crc65$CCLE[, ii, drop = FALSE], sites = m.crc65[, ii])
ll$medico <- runMod(drugMat = drug.crc65$MEDICO[, ii, drop = FALSE], sites = m.crc65[, ii])
ll$gdsc <- runMod(drugMat = drug.crc65$GDSC[, ii, drop = FALSE], sites = m.crc65[, ii])
ll$ctrp <- runMod(drugMat = drug.crc65$CTRP[, ii, drop = FALSE], sites = m.crc65[, ii])

saveRDS(ll, file = "Res/20170731_weight.comboScore/elnet.allModules.crc65_alldrugs.RDS")