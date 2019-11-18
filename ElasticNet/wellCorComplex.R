library(ggplot2)
library(openxlsx)
site.tryp <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170110_mq15processLOD2/sites.trypsin.RDS")
site.gluc <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170110_mq15processLOD2/sites.gluC.RDS")
prot.tryp.nci60 <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
prot.tryp.crc65 <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170110_mq15processLOD2/crc65.proteinGroups.RDS")
prot.gluc.nci60 <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170719_nci60.proteinGroups.gluc/nci60.proteinGroups.gluc.RDS")

res <- readRDS(file = "Res/20170926_corumCorrelation/summarizedRes.RDS")

annotPlot <- function(r, expr, annot, xref, complex, alpha = 0.01) {
  subr <- r[r$complex == complex, ]
  i <- which(subr$pval < 0.01)
  if (length(i) < 1)
    return(NULL)
  
  
  wc1 <- unlist(strsplit(subr$feature1[i], ";"))
  wc2 <- unlist(strsplit(subr$feature2[i], ";"))
  
  af <- unique(c(wc1, wc2))
  an <- annot[af, c("Protein", "Protein names", "Gene names", "Position", "Localization prob", 
                    "Regulatory site", "Regulatory site function", "Sequence window")]
  
  ft <- rep(paste(wc1, wc2, sep = "_"), each = ncol(expr))
  ft <- gsub(" ", ".", ft)
  if (is.null(xref$too.short))
    xref$too.short <- "CO"
  df <- data.frame(psite1 = c(t(expr[wc1, ])),
                   psite2 = c(t(expr[wc2, ])),
                   fac = ft, 
                   col = rep(xref$too.short, length(wc1)))
  
  df <- df[!is.na(df$psite1)& !is.na(df$psite2), ]
  p <- ggplot(data = df, aes(x = psite2, y = psite1, color = col)) + geom_point()
  p <- p + facet_wrap(~fac)
  p <- p + geom_abline(slope=1, intercept=0)
  
  list(annot = an, plot = p)
}


aplot <- function(complex, alpha, data) {
  if (data == "nci60 trypsin") {
    r = res$site.trypsin.nci60
    expr = site.tryp$nci60_expr
    annot = site.tryp$annot
    xref = site.tryp$nci60_xref
  } else if (data == "nci60 gluc") {
    r = res$site.gluc.nci60
    expr = site.gluc$nci60_expr
    annot =  site.gluc$annot
    xref = site.gluc$nci60_xref
  } else if (data == "crc65") {
    cl <- setdiff(colnames(site.tryp$crc65_expr), "CoCM-1_CRC65")
    r = res$site.trypsin.crc65
    expr = site.tryp$crc65_expr[, cl]
    annot = site.tryp$annot
    xref = site.tryp$crc65_xref[cl, ]
  } else {
    stop("data should be one of 'nci60 trypsin','nci60 gluc' or 'crc65'.")
  }
  annotPlot(complex = complex, alpha = alpha, r = r, expr = expr, annot = annot, xref = xref)
}

## ======================================= NCI60 TRYPSIN ==========================================
## MAPK3, MAPK1
r <- aplot(alpha = 0.01, complex = "VEGFR2-S1PR5-ERK1/2-PKC-alpha complex", data = "nci60 trypsin")
write.xlsx(r$annot, file = "Res/20170927_corumCorrelationPlots/wellcor_nci60.trypsin_MAPK3.MAPK1.xlsx", asTable = TRUE)
png("Res/20170927_corumCorrelationPlots/wellcor_nci60.trypsin_MAPK3.MAPK1.png", width = 10, height = 8, units = "cm", res = 300)
r$plot
dev.off()

# HDAC1, HDAC2
r <- aplot(alpha = 0.01, complex = "DNTTIP1-ZNF541-HDAC1-HDAC2 complex", data = "nci60 trypsin")
r <- aplot(alpha = 0.01, complex = "BRCA1-HDAC1-HDAC2 complex", data = "nci60 trypsin")
r <- aplot(alpha = 0.01, complex = "MiDAC complex", data = "nci60 trypsin")
r$annot[1]
write.xlsx(r$annot, file = "Res/20170927_corumCorrelationPlots/wellcor_nci60.trypsin_HDAC1.HDAC2.xlsx", asTable = TRUE)
png("Res/20170927_corumCorrelationPlots/wellcor_nci60.trypsin_HDAC1.HDAC2.png", width = 10, height = 8, units = "cm", res = 300)
r$plot
dev.off()


# ORC1, CDC6
r <- aplot(alpha = 0.01, complex = "BRCA1-IRIS-pre-replication complex", data = "nci60 trypsin")
r$annot[1]
write.xlsx(r$annot, file = "Res/20170927_corumCorrelationPlots/wellcor_nci60.trypsin_CDC6.ORC1.xlsx", asTable = TRUE)
png("Res/20170927_corumCorrelationPlots/wellcor_nci60.trypsin_CDC6.ORC1.png", width = 18, height = 14, units = "cm", res = 300)
r$plot
dev.off()


# CDK2, TFDP1
r <- aplot(alpha = 0.01, complex = "E2F-1-DP-1-cyclinA-CDK2 complex", data = "nci60 trypsin")
r$annot[1]
write.xlsx(r$annot, file = "Res/20170927_corumCorrelationPlots/wellcor_nci60.trypsin_CDK2.TFDP1.xlsx", asTable = TRUE)
png("Res/20170927_corumCorrelationPlots/wellcor_nci60.trypsin_CDK2.TFDP1.png", width = 10, height = 8, units = "cm", res = 300)
r$plot
dev.off()


## ======================================= CRC65 TRYPSIN ==========================================
## MAPK3, MAPK1
r <- aplot(alpha = 0.01, complex = "VEGFR2-S1PR5-ERK1/2-PKC-alpha complex", data = "crc65")
r$annot[1]
write.xlsx(r$annot, file = "Res/20170927_corumCorrelationPlots/wellcor_crc65.trypsin_MAPK3.MAPK1.xlsx", asTable = TRUE)
png("Res/20170927_corumCorrelationPlots/wellcor_crc65.trypsin_MAPK3.MAPK1.png", width = 10, height = 8, units = "cm", res = 300)
r$plot
dev.off()


# YEATS ... histone acetylation
r <- aplot(alpha = 0.01, complex = "ATAC complex, GCN5-linked", data = "crc65")
r <- aplot(alpha = 0.01, complex = "ATAC complex, YEATS2-linked", data = "crc65")
r$annot[1]
write.xlsx(r$annot, file = "Res/20170927_corumCorrelationPlots/wellcor_crc65.trypsin_YEATS.xlsx", asTable = TRUE)
png("Res/20170927_corumCorrelationPlots/wellcor_crc65.trypsin_YEATS.png", width = 18, height = 14, units = "cm", res = 300)
r$plot
dev.off()


## ======================================= NCI60 GLUC ==========================================
r <- aplot(alpha = 0.01, complex = "Polycomb repressive complex 1 (PRC1, hPRC-H)", data = "nci60 gluc")
r$annot[1]
r$plot

# BRCA1, BARD1
r <- aplot(alpha = 0.01, complex = "BRCA1-BARD1-UbcH7c complex", data = "nci60 gluc")
r <- aplot(alpha = 0.01, complex = "BRCA1 C complex", data = "nci60 gluc")
r$annot[1]
write.xlsx(r$annot, file = "Res/20170927_corumCorrelationPlots/wellcor_nci60.gluc_BRCA1.BARD1.xlsx", asTable = TRUE)
png("Res/20170927_corumCorrelationPlots/wellcor_nci60.gluc_BRCA1.BARD1.png", width = 18, height = 14, units = "cm", res = 300)
r$plot
dev.off()

## ================================== 
source("R/drugResponseProcess/extract.drugmat.nci60.R")
r <- aplot(alpha = 0.01, complex = "DNTTIP1-ZNF541-HDAC1-HDAC2 complex", data = "nci60 trypsin")
emat <- site.tryp$nci60_expr[rownames(r$annot), ]
dmat <- drug.nci60$dtp.gi50
cc <- cor(t(emat), t(dmat), use = "pair")

which(t(cc) > 0.35, arr.ind = TRUE)
hist(cc, breaks = 100)





## 

p.expr <- site.tryp$crc65_expr
p.annot <- site.tryp$annot
p.xref <- sites.try 

g1 <- "53283 CDK1"
g2 <- "53960 CDK2"
an <- p.annot[c(g1, g2), c("Protein", "Protein names", "Gene names", "Position", "Localization prob", 
                     "Regulatory site", "Regulatory site function", "Sequence window")]
write.xlsx(an, file = "Res/20170927_corumCorrelationPlots/wellcor_crc65.trypsin_CDK1.CDK2.xlsx", asTable = TRUE)
png("Res/20170927_corumCorrelationPlots/wellcor_crc65.trypsin_CDK1.CDK2.png", width = 12, height = 12, units = "cm", res = 300)
plot(p.expr[g1, ], p.expr[g2, ], xlab = g1, ylab = g2, col = "orange2", pch = 19)
abline(a =0, b = 1)
dev.off()

