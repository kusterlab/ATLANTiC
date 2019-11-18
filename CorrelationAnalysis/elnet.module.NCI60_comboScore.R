library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)
library(psych)
library(reshape2)

source("../R/elnetFit.R")

# drug combination scores
dat <- readRDS("Res/20170731_weight.comboScore/weighted.scores.RDS")
d.mat <- dat$weighted.score[dat$filter, ]
h <- dat$comboScore[dat$filter, 1:6]
rownames(d.mat) <- paste(h$NAME1, h$CONC1, h$NAME2, h$CONC2, sep = "_")
colnames(d.mat) <- gsub("-| |\\(TB\\)|/|ATCC", "", colnames(d.mat))
colnames(d.mat)[colnames(d.mat) == "UO31"] <- "U031"
colnames(d.mat)[colnames(d.mat) == "NCIADRRES"] <- "NCIADRES"
colnames(d.mat)[colnames(d.mat) == "7860"] <- "786O"
colnames(d.mat) <- paste0(colnames(d.mat), "_NCI60")

### input mat
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


setdiff(colnames(d.mat), colnames(m.nci60))
setdiff(colnames(m.nci60), colnames(d.mat))




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
res <- runMod(drugMat = d.mat[, colnames(m.nci60)], sites = m.nci60)
saveRDS(res, file = "Res/20170731_weight.comboScore/elnet.allModules.nci60_drugCombo.RDS")
