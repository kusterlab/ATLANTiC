setwd("/media/general/projects/NCI_60_phospho/Chen/")

library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)
library(omic3plus)
library(impute)
library(WebGestaltR)

dat <- readRDS("Res/20170904_concordanceData/nci60_dtp.RDS")
dproc <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")

icell <- dat$xref$too.short != "LE"
y <- dat$y$data[, icell]
X <- lapply(dat$x, function(x) x$imputed[, icell])

y <- y[rowMaxs(y) > 1.5, ]
X <- lapply(X, function(x) {
  x[rowSds(x)>1e-5, ]
  # x[rowSds(x)>1, ]
})

sapply(X, dim)

cho.ky <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
cho.kx <- 50*(2^(0:3))

vp <- concord(X, y, kx = cho.kx, ky = cho.ky, option = "nk", ncomp = 20,
              scan = FALSE, ncores = 5,
              dmod = 1, center.x = TRUE, scale.x = TRUE,
              center.y = TRUE, scale.y = TRUE, pos = FALSE)

loadingy <- data.frame(vp$loading.y)
loadingy <- cbind(drug = rownames(y), dat$drugannot[rownames(y), "MOA"], loadingy)

saveRDS(list(input = dat, res = vp), file = "Res/20171107_concord2/nci60_dtp_ky10.RDS")


rownames(vp$loading.x)[vp$loading.x[, 26] < 0]

#### PC - drug association
# PC01 neg - DNA(T2) inhibitors - (R) ABCB1
# PC02 neg - DNA(A7) inhibitors - (S) (TSNAX, PFKM)
# PC03 neg - DNA inhibitors -
# PC04 neg - multi-kinase inhibitors - (S) MYLK, FLNC, TPM1, PRKCA, PPP1R12A
# PC07 neg - VEGFR, FGFR hibibitors (R) (HDGFRP2, PAICS) (S) (DCLK2)
# PC09 neg - DNA inihibtors - ribosome enriched
# PC10 neg - HMG-CoA reductase inhibitor - PRKAA1
# PC11 neg - Hormone (Estrogen) inhibitor - (S) GREB1*,NCOA3(AIB1),RARA,RBM39,BCAS3
# PC12 neg - DNA inhibitors
# PC15 neg - DNA inhibitors - (R) TP53BP1 (pathway: DNA Double-Strand Break Repair and ATM Pathway)
# PC16 pos - STK,YK - ALK inhibitors - NTRK1, SOS1, SH2B1
# PC17 pos - STK, YK - (S) STK4
# PC19 neg - DNA inhibitors (anti-metabolites)
# PC20 neg - DNA(A7) inhibitors
# PC21 neg - hormone inhibitor (AR) - (S) CALU, NUCB1, CALD1 (ER/calcium homeostasis), PELP1 (ER)
# PC23 neg - EGFR inhibitors - GAB1, GRB7, PRKCI, PTPN11
# PC26 neg - Db (S) TPD52L2

dfl <- lapply(1:30, function(i) {
  df <- data.frame (id = rownames(vp$loading.x)[vp$loading.x[, i] != 0], 
                    val = vp$loading.x[vp$loading.x[, i] != 0, i])
  df <- df[order(df$val), ]
  df$PC <- i
  df
})
dflm <- do.call(rbind, dfl)

i <- 3
aa <- unique(sapply(strsplit(as.character(dfl[[i]]$id), " |;"), "[", 2))
write.table(data.frame(aa), file = "tmp.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
tt <- dfl[[i]]


# library(WebGestaltR)
# all <- unique(na.omit(sapply(strsplit(rownames(vp$loading.x), " |;"), "[", 2)))
# ll <- lapply(1:30, function(i) {
#   itx <- unique(na.omit(sapply(strsplit(rownames(vp$loading.x)[vp$loading.x[, i] != 0], " |;"), "[", 2)))  
#   WebGestaltR(enrichMethod="ORA", organism="hsapiens", enrichDatabase = "pathway_Wikipathway",
#               interestGene = itx, referenceGene = all, interestGeneType = "genesymbol",
#                       referenceGeneType = "genesymbol", is.output = FALSE, fdrThr = 1)
# })
# 
# nr <- sapply(sapply(ll, nrow), function(x) ifelse(is.null(x), 0, x))
# names(nr) <- paste0("PC", 1:30)
# nr
# tt <- ll[[26]]
# plot(vp$score.x[, 11], vp$score.y[, 11])




