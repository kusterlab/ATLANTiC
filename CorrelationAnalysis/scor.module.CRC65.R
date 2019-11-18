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

ii <- intersect(colnames(drug.crc65$CCLE), colnames(m.crc65))

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
ll$ccle <- scor(dmat = drug.crc65$CCLE[, ii, drop = FALSE], pcmat = m.crc65[, ii])
ll$medico <- scor(dmat = drug.crc65$MEDICO[, ii, drop = FALSE], pcmat = m.crc65[, ii])
ll$gdsc <- scor(dmat = drug.crc65$GDSC[, ii, drop = FALSE], pcmat = m.crc65[, ii])
ll$ctrp <- scor(dmat = drug.crc65$CTRP[, ii, drop = FALSE], pcmat = m.crc65[, ii])

saveRDS(ll, file = "Res/20170731_weight.comboScore/scor.allModules.crc65_alldrugs.RDS")




