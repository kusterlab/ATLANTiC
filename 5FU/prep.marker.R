library(randomForest)
library(rpart)
library(matrixStats)

source("/media/general/projects/NCI_60_phospho/Chen/R/drugResponseProcess/extract.drugmat.nci60.R")
pgs <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
pgs.crc <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170110_mq15processLOD2/crc65.proteinGroups.RDS")
tryp.phos <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170110_mq15processLOD2/sites.trypsin.RDS")
gluc.phos <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170110_mq15processLOD2/sites.gluC.RDS")

## drug sensitivity data
grep("Fluorouracil", drug.nci60$dtp.annot$Name)
fu19893 <- drug.nci60$dtp.gi50[20, ]
fu757036 <- drug.nci60$dtp.gi50[214, ]
plot(fu19893, fu757036)

markers <- c("TYMP", "TK1", # 5fu -> FUDR -> FdUMP
             "UMPS", "PPAT", # 5fu -> FUMP 
             "UPP1",  # 5fu -> FUR, "UPP2" as well, but not detected
             "UCK2", "UCK1", # FUR -> FUMP
             "RRM1", "RRM2", # FUDP ->FdUTP
             "SLC29A1", # importer, "SLC22A7" as well, but not detected
             "ABCC3", "ABCC4", # exporter, "ABCC5", "ABCG2" as well, but not detected
             "TYMS", # target
             "DPYD") # 5fu -> DHFU (degradation)

dat <- list()
dat$fu19893 <- fu19893
dat$fu757036 <- fu757036

dat$crc.xref <- pgs.crc$xref
dat$nci.xref <- pgs$xref
## NCI 60
ir <- pgs$annot$`Gene names` %in% markers
expr <- pgs$ibaq[ir, ]
expr <- t(expr)
colnames(expr) <- sapply(strsplit(colnames(expr), " "), "[", 2)
expr <- data.frame(expr)
dat$nci60.prot <- expr

## CRC65
ir <- pgs.crc$annot$`Gene names` %in% markers
expr <- pgs.crc$ibaq[ir, ]
expr <- t(expr)
colnames(expr) <- sapply(strsplit(colnames(expr), " "), "[", 2)
expr <- data.frame(expr)
dat$crc65.prot <- expr


## NCI 60 phospho TK1
ir <- which(tryp.phos$annot$`Gene names` == "TK1")
an <- tryp.phos$annot[ir, ]
expr <- data.frame(t(tryp.phos$nci60_expr[ir, ]))
colnames(expr) <- paste(an$`Gene name`, substr(an$`Sequence window`, 16, 16), 
                        an$`Positions within proteins`, sep = "_")
dat$nci60.tk1phos <- expr


## CRC 65 phospho TK1
expr <- data.frame(t(tryp.phos$crc65_expr[ir, ]))
colnames(expr) <- paste(an$`Gene name`, substr(an$`Sequence window`, 16, 16), 
                        an$`Positions within proteins`, sep = "_")
dat$crc65.tk1phos <- expr

## NCI 60 phos TK1, gluc
ir <- which(gluc.phos$annot$`Gene name` == "TK1")
an <- gluc.phos$annot[ir, ]
expr <- data.frame(v =gluc.phos$nci60_expr[ir, ])
colnames(expr) <- paste(an$`Gene name`, substr(an$`Sequence window`, 16, 16), 
                        an$`Positions within proteins`, sep = "_")
dat$nci60.tk1phos.gluc <- expr

## binarize data
dat$bin <- function(x1, x2, cut) {
  plot(x1, x2)
  abline(v = cut, h = cut)
  xx <- as.numeric(x1 > cut & x2 > cut)
  xx[is.na(xx)] <- 0
  xx
}


## imputed
dat$lod2 <- function(x) {
  x[seq_along(x)] <- lapply(x, function(xs) {
    xs[is.na(xs)] <- min(xs, na.rm = TRUE) - log10(2)
    xs
  })
  x
}

summary(dat)
saveRDS(dat, file = "Res/20170317_Fluorouracil/prep.dat.RDS")
