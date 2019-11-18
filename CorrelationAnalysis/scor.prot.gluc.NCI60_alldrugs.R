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
scor <- function(dmat, pcmat) {
  
  cat("Calculating correlation ... \n")
  mcc <- mclapply(1:nrow(dmat), function(i) {
    print(i)
    res <- corr.test(t(pcmat), t(dmat[i, , drop = FALSE]), adjust = "none", use = "pair")
    df <- melt(res$r)
    df$n <- c(res$n)
    df$t <- c(res$t)
    df$p <- c(res$p)
    df$se <- c(res$se)
    df <- cbind(df, res$ci)
    return(df)
  }, mc.cores = 10)
  # cat("Filtering resulsts ... \n")
  # mcc <- lapply(mcc, function(x) x[ abs(x$value) > 0.25, ])
  cat("Binding resulsts ... \n")
  do.call(rbind, mcc)
  
}

expr <- x$fun_prepElnet(x, celllines = "all")
dim(expr$imputed)
exprm <- expr$imputed
colnames(exprm) <- gsub("_GluC", "_NCI60", colnames(exprm))
identical(colnames(exprm), colnames(drug.nci60$ccle.auc))

ll <- list()
ll$ccle <- scor(dmat = drug.nci60$ccle.auc, pcmat = exprm)
ll$dtp <- scor(dmat = drug.nci60$dtp.gi50, pcmat = exprm)
ll$gdsc <- scor(dmat = drug.nci60$gdsc.auc, pcmat = exprm)
ll$ctrp <- scor(dmat = drug.nci60$ctrp.auc, pcmat = exprm)

saveRDS(ll, file = "Res/20170731_weight.comboScore/scor.prot.gluc.nci60_imputed_alldrugs.RDS")

### no imputation


exprm_nf <- x$ibaq[rownames(exprm), ]
colnames(exprm_nf) <- colnames(exprm)

ll <- list()
ll$ccle <- scor(dmat = drug.nci60$ccle.auc, pcmat = exprm_nf)
ll$dtp <- scor(dmat = drug.nci60$dtp.gi50, pcmat = exprm_nf)
ll$gdsc <- scor(dmat = drug.nci60$gdsc.auc, pcmat = exprm_nf)
ll$ctrp <- scor(dmat = drug.nci60$ctrp.auc, pcmat = exprm_nf)

saveRDS(ll, file = "Res/20170731_weight.comboScore/scor.prot.gluc.nci60_noimputation_alldrugs.RDS")
