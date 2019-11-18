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
x <- readRDS(file = "Res/20170719_nci60.proteinGroups.gluc/nci60.proteinGroups.gluc.RDS")


###
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


expr <- x$fun_prepElnet(x, celllines = "all")
dim(expr$imputed)
exprm <- expr$imputed
colnames(exprm) <- gsub("_GluC", "_NCI60", colnames(exprm))
identical(colnames(exprm), colnames(drug.nci60$ccle.auc))

## all phospho-sites
ll <- list()
ll$ccle <- runMod(drugMat = drug.nci60$ccle.auc, sites = exprm)
ll$dtp <- runMod(drugMat = drug.nci60$dtp.gi50, sites = exprm)
ll$gdsc <- runMod(drugMat = drug.nci60$gdsc.auc, sites = exprm)
ll$ctrp <- runMod(drugMat = drug.nci60$ctrp.auc, sites = exprm)

saveRDS(ll, file = "Res/20170731_weight.comboScore/elnet.prot.gluc.nci60_alldrugs.RDS")
