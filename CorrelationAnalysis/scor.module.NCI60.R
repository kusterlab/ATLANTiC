library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)
library(psych)
library(reshape2)

setwd("/media/general/projects/NCI_60_phospho/Chen")
library(reshape2)
library(matrixStats)
library(parallel)
library(psych)

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

identical(colnames(m.nci60), colnames(drug.nci60$ctrp.auc))

###
scor <- function(dmat, pcmat) {
  
  cat("Calculating correlation ... \n")
  mcc <- mclapply(1:nrow(dmat), function(i) {
    print(i)
    res <- corr.test(t(pcmat), t(dmat[i, , drop = FALSE]), adjust = "none")
    df <- melt(res$r)
    df$n <- c(res$n)
    df$t <- c(res$t)
    df$p <- c(res$p)
    df$se <- c(res$se)
    df <- cbind(df, res$ci)
    return(df)
  }, mc.cores = 10)
  cat("Filtering resulsts ... \n")
  mcc <- lapply(mcc, function(x) x[ abs(x$value) > 0.25, ])
  cat("Binding resulsts ... \n")
  do.call(rbind, mcc)
  
}

ll <- list()
ll$dtp <- scor(drug.nci60$dtp.gi50, m.nci60)
ll$ccle <- scor(drug.nci60$ccle.auc, m.nci60)
ll$gdsc <- scor(drug.nci60$gdsc.auc, m.nci60)
ll$ctrp <- scor(drug.nci60$ctrp.auc, m.nci60)

saveRDS(ll, file = "Res/20170731_weight.comboScore/scor.allModules.nci60_alldrugs.RDS")

















