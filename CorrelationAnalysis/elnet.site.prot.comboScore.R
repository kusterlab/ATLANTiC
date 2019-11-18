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
site.tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
site.tryp.mat <- site.tryp$fun_prepElnet(site.tryp, panel = "NCI60", celllines = "all")

site.gluc <- readRDS("Res/20170110_mq15processLOD2/sites.gluC.RDS")
site.gluc.mat <- site.gluc$fun_prepElnet(site.gluc, celllines = "all")

# prot.tryp <- readRDS("Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
# prot.tryp.mat <- prot.tryp$fun_prepElnet(prot.tryp, celllines = "all")
# 
# prot.gluc <- readRDS("Res/20170719_nci60.proteinGroups.gluc/nci60.proteinGroups.gluc.RDS")
# prot.gluc.mat <- prot.gluc$fun_prepElnet(prot.gluc, celllines = "all")


identical(colnames(d.mat), colnames(site.tryp.mat$imputed))
all(colnames(d.mat) %in% colnames(site.tryp.mat$imputed))

## Define function

runMod <- function( drugMat, sites, regsitesOnly = FALSE, prefix = "") {# the drug matrix, rows are drugs columns are samples
  
  if (! identical(colnames(drugMat), colnames(sites)))
    stop("check columns order.")
  
  fit <- list()
  for (i in 1:nrow(drugMat)) {
    if (i %% 100 == 0)
      gc()
    print(paste(prefix, i))
    gi50 <- drugMat[i, ]
    iSample <- !is.na(gi50)
    
    if (sum(iSample) < 12) {
      fit[[i]] <- NA
      next()
    }
    
    xmat <- sites[, iSample]
    fmod <- try(elnetFit(t(xmat), 
                             y = gi50[iSample], 
                             standardizex = TRUE,
                             nboot = 100, 
                             permpvalue = 0, 
                             x.na.index = matrix(FALSE, nrow(xmat), ncol(xmat)),
                             times=1, 
                             cores = 16, 
                             write.table = FALSE), silent = TRUE)
    if (inherits(fmod, "list")) {
      ss <- fmod$summary
      fit[[i]] <- ss[which(ss$nonzerocounts > 5), ]
    } else
      fit[[i]] <- NA
  }
  names(fit) <- rownames(drugMat)
  fit
}


# x <- prot.tryp.mat$imputed
# res <- runMod(drugMat = d.mat[, colnames(x)], sites = x, prefix = "prot.tryp")
# saveRDS(res, file = "Res/20170810_comboScoreElnet/elnet.prot.tryp.nci60_drugCombo.RDS")
# rm(res)
# gc()
# 
# x <- prot.gluc.mat$imputed
# colnames(x) <- gsub("_GluC", "_NCI60", colnames(x))
# res <- runMod(drugMat = d.mat[, colnames(x)], sites = x, prefix = "prot.gluc")
# saveRDS(res, file = "Res/20170810_comboScoreElnet/elnet.prot.gluc.nci60_drugCombo.RDS")
# rm(res)
# gc()

#
# x <- site.tryp.mat$imputed
# res <- runMod(drugMat = d.mat[, colnames(x)], sites = x, prefix = "site.tryp")
# saveRDS(res, file = "Res/20170810_comboScoreElnet/elnet.site.tryp.nci60_drugCombo.RDS")
# rm(res)
# gc()

###
x <- site.gluc.mat$imputed
colnames(x) <- paste0(colnames(x), "_NCI60")
res <- runMod(drugMat = d.mat[, colnames(x)], sites = x, prefix = "site.gluc")
saveRDS(res, file = "Res/20170810_comboScoreElnet/elnet.site.gluc.nci60_drugCombo.RDS")
rm(res)
gc()

