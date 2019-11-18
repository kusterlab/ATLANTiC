library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)
library(psych)
library(reshape2)

# source("../R/elnetFit.R")

# ====================== drug combination scores =========================
dat <- readRDS("Res/20170731_weight.comboScore/weighted.scores.RDS")
d.mat <- dat$weighted.score[dat$filter, ]
h <- dat$comboScore[dat$filter, 1:6]
rownames(d.mat) <- paste(h$NAME1, h$CONC1, h$NAME2, h$CONC2, sep = "_")
colnames(d.mat) <- gsub("-| |\\(TB\\)|/|ATCC", "", colnames(d.mat))
colnames(d.mat)[colnames(d.mat) == "UO31"] <- "U031"
colnames(d.mat)[colnames(d.mat) == "NCIADRRES"] <- "NCIADRES"
colnames(d.mat)[colnames(d.mat) == "7860"] <- "786O"
colnames(d.mat) <- paste0(colnames(d.mat), "_NCI60")

# ====================== load molecular features =========================

site.tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
site.tryp.mat <- site.tryp$fun_prepElnet(site.tryp, panel = "NCI60", celllines = "all")

site.gluc <- readRDS("Res/20170110_mq15processLOD2/sites.gluC.RDS")
site.gluc.mat <- site.gluc$fun_prepElnet(site.gluc, celllines = "all")

prot.tryp <- readRDS("Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
prot.tryp.mat <- prot.tryp$fun_prepElnet(prot.tryp, celllines = "all")

prot.gluc <- readRDS("Res/20170719_nci60.proteinGroups.gluc/nci60.proteinGroups.gluc.RDS")
prot.gluc.mat <- prot.gluc$fun_prepElnet(prot.gluc, celllines = "all")


# validate colnames
identical(colnames(site.tryp.mat$imputed), paste0(colnames(site.gluc.mat$imputed), "_NCI60"))
identical(colnames(site.tryp.mat$imputed), colnames(prot.tryp.mat$imputed))
identical(colnames(site.tryp.mat$imputed), gsub("_GluC", "_NCI60", colnames(prot.gluc.mat$imputed)))

# reorder d.mat columns
d.mat <- d.mat[, colnames(site.tryp.mat$imputed)]
identical(colnames(site.tryp.mat$imputed), colnames(d.mat))

# ========================== calculation =====================================

rmad <- rowMads(d.mat, na.rm = TRUE)
names(rmad) <- rownames(d.mat)

scorr <- function(x, dmat, vmad) {
  mcc <- mclapply(1:nrow(dmat), function(i) {
    if (i %% 100 == 0 || i > 4400)
      print(i)
    res <- corr.test(t(x), t(dmat[i, , drop = FALSE]), adjust = "none", use = "pair")
    df <- melt(res$r)
    df$n <- c(res$n)
    df$t <- c(res$t)
    df$p <- c(res$p)
    df$se <- c(res$se)
    df <- cbind(df, res$ci)
    return(df)
  }, mc.cores = 28)
  
  vv <<- mcc
  mccf <- lapply(mcc, function(x) {
    if (is.null(x))
      return(NULL)
    x[ abs(x$value) > 0.3, ]
    })
  mcc1 <- do.call(rbind, mccf)
  mcc1$scoreMad <- vmad[mcc1$Var2]
  mcc1
}


# cat("site tryp \n")
# cat(as.character(Sys.time()))
# ll <- scorr(site.tryp.mat$imputed, d.mat, vmad = rmad)
# saveRDS(ll, file = "Res/20170801_comboScoreFittings/scor.imputed.comboScore_site.tryp.RDS")
# rm(ll)
# gc()

# cat("site gluc \n")
# cat(as.character(Sys.time()))
# ll <- scorr(site.gluc.mat$imputed, d.mat, vmad = rmad)
# saveRDS(ll, file = "Res/20170801_comboScoreFittings/scor.imputed.comboScore_site.gluc.RDS")
# rm(ll)
# gc()

# 
# cat("prot tryp \n")
# cat(as.character(Sys.time()))
# ll <- scorr(prot.tryp.mat$imputed, d.mat, vmad = rmad)
# saveRDS(ll, file = "Res/20170801_comboScoreFittings/scor.imputed.comboScore_prot.tryp.RDS")
# rm(ll)
# gc()


# cat("prot gluc \n")
# cat(as.character(Sys.time()))
# ll <- scorr(prot.gluc.mat$imputed, d.mat, vmad = rmad)
# saveRDS(ll, file = "Res/20170801_comboScoreFittings/scor.imputed.comboScore_prot.gluc.RDS")
# rm(ll)
# gc()
# 
# 
# 
# #######
# 
removeimp <- function(x, i) {
  x[i] <- NA
  x
}
# 
# cat("site tryp \n")
# cat(as.character(Sys.time()))
# ll <- scorr(removeimp(site.tryp.mat$imputed, site.tryp.mat$na_index), d.mat, vmad = rmad)
# saveRDS(ll, file = "Res/20170801_comboScoreFittings/scor.noimputation.comboScore_site.tryp.RDS")
# rm(ll)
# gc()
# 
# 
# cat("site gluc \n")
# cat(as.character(Sys.time()))
# ll <- scorr(removeimp(site.gluc.mat$imputed, site.gluc.mat$na_index), d.mat, vmad = rmad)
# saveRDS(ll, file = "Res/20170801_comboScoreFittings/scor.noimputation.comboScore_site.gluc.RDS")
# rm(ll)
# gc()


cat("prot tryp \n")
cat(as.character(Sys.time()))
ll <- scorr(x = removeimp(prot.tryp.mat$imputed, prot.tryp.mat$na_index), d.mat, vmad = rmad)
saveRDS(ll, file = "Res/20170801_comboScoreFittings/scor.noimputation.comboScore_prot.tryp.RDS")
rm(ll)
gc()


cat("prot gluc \n")
cat(as.character(Sys.time()))
ll <- scorr(removeimp(prot.gluc.mat$imputed, prot.gluc.mat$na_index), d.mat, vmad = rmad)
saveRDS(ll, file = "Res/20170801_comboScoreFittings/scor.noimputation.comboScore_prot.gluc.RDS")
rm(ll)
gc()















