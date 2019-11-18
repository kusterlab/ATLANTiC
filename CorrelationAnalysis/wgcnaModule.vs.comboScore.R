library(psych)
library(parallel)
library(reshape2)
library(matrixStats)
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


####
dat <- readRDS("Res/20170731_weight.comboScore/weighted.scores.RDS")
d.mat <- dat$weighted.score[dat$filter, ]
h <- dat$comboScore[dat$filter, 1:6]
rownames(d.mat) <- paste(h$NAME1, h$CONC1, h$NAME2, h$CONC2, sep = "_")
colnames(d.mat) <- gsub("-| |\\(TB\\)|/|ATCC", "", colnames(d.mat))
colnames(d.mat)[colnames(d.mat) == "UO31"] <- "U031"
colnames(d.mat)[colnames(d.mat) == "NCIADRRES"] <- "NCIADRES"
colnames(d.mat)[colnames(d.mat) == "7860"] <- "786O"
colnames(d.mat) <- paste0(colnames(d.mat), "_NCI60")

setdiff(colnames(d.mat), colnames(m.nci60))
setdiff(colnames(m.nci60), colnames(d.mat))

rmad <- rowMads(d.mat, na.rm = TRUE)
names(rmad) <- rownames(d.mat)

########
cids <- colnames(m.nci60)

mcc <- mclapply(1:nrow(d.mat), function(i) {
  print(i)
  res <- corr.test(t(m.nci60[, cids]), t(d.mat[i, cids, drop = FALSE]), adjust = "none")
  df <- melt(res$r)
  df$n <- c(res$n)
  df$t <- c(res$t)
  df$p <- c(res$p)
  df$se <- c(res$se)
  df <- cbind(df, res$ci)
  return(df)
}, mc.cores = 20)

mccf <- lapply(mcc, function(x) x[ abs(x$value) > 0.3, ])
mcc1 <- do.call(rbind, mccf)
mcc1$scoreMad <- rmad[mcc1$Var2]

saveRDS(mcc1, file = "Res/20170731_weight.comboScore/module.drugCombo.RDS")


length(unique(mcc1$scoreMad))
length(unique(mcc1$Var2))








