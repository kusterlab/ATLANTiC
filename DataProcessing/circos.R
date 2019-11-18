library(circlize)
source("R/final/circos/circos.trackBarplot.R")

### ===========================================================================================
###
###                         NCI60 and CRC65 
###
### ===========================================================================================
tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
crc.phos <- readRDS("Res/20170721_numbersFinalPlotData/numberid_crc65_PHOS.RDS")
crc.fp <- readRDS("Res/20170721_numbersFinalPlotData/numberid_crc65_FP.RDS")
nci.phos <- readRDS("Res/20170721_numbersFinalPlotData/numberid_nci60_PHOS.RDS")
nci.fp <- readRDS("Res/20170721_numbersFinalPlotData/numberid_nci60_FP.RDS")

transparentColor <- function(x, alpha = 100) {
  rgbs <- col2rgb(x)
  rgb(rgbs[1, ], rgbs[2, ], rgbs[3, ], alpha = 100, maxColorValue = 255)
}

# 
xref <- tryp$nci60_xref
xref <- xref[order(xref$too.short), ]
ccol <- unique(xref$color[xref$too.short == "CO"])

##  fullprotein
cellline <- cbind(nci.fp$cellline, rbind(0, crc.fp$cellline, 0))
wid <- c(nci.fp$tissue["nsample", ], crc.fp$tissue["nsample", ])
fcc <- c(unique(xref$color), ccol)
cc <- transparentColor(fcc)
tissue <- c(nci.fp$tissue["tisid", ], CRC65 = unname(crc.fp$tissue["tisid", ]))

tot <- rep(tissue, wid)
fp <- rbind(cellline, tot - colSums(cellline))
      
## phospho-site
cellline <- cbind(nci.phos$cellline, rbind(0, crc.phos$cellline, 0))
wid <- c(nci.phos$tissue["nsample", ], crc.phos$tissue["nsample", ])
fcc <- c(unique(xref$color), ccol)
cc <- transparentColor(fcc)
tissue <- c(nci.phos$tissue["tisid", ], CRC65 = unname(crc.phos$tissue["tisid", ]))

tot <- rep(tissue, wid)
phos <- rbind(cellline, tot - colSums(cellline))

fp <- apply(fp, 2, cumsum)
phos <- apply(phos, 2, cumsum)

c1 <- rbind(rep("gray25", ncol(phos)), 
            rep("gray50", ncol(phos)), 
            rep("gray75", ncol(phos)), 
            rep(cc, wid))

numlist <- list(protein = fp, phos = phos)
colorlist <- list(c1, c1)

### ===========================================================================================
###
###                         NCI60 and CRC65 
###
### ===========================================================================================
ss <- readRDS(file = "Res/20170307_countingComparison/cnt.OV.BR.hela.nci60.crc65.RDS")
sscol <- list(c("blue", "gray"), c("blue", "gray"), c("blue"), 
            rbind(rep("blue", 60), rep("orange", 60), xref$color, rep("gray", 60)), 
            c("blue", "gray"))
names(sscol) <- names(ss)

sscol <- sscol[1:3]
ss <- ss[1:3]




### ===========================================================================================
###
###                         NCI60 and CRC65 
###
### ===========================================================================================


# everything together
l <- c(ss, numlist)
c <- c(sscol, colorlist)
pdf(file = "Res/20170721_numbersFinalPlotData/circo_all.pdf", width = 7, height = 7)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.barplot(l = l, col = c, bg.border = "white", track.height = 0.6)
dev.off()


# only phosphosite
l <- c(ss, numlist[-1])
c <- c(sscol, colorlist[-1])
pdf(file = "Res/20170721_numbersFinalPlotData/circo_all_psite_includeExternal.pdf", width = 7, height = 7)
# circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.barplot(l = l, col = c, bg.border = "white", track.height = 0.6)
dev.off()

# only phosphosite this study
l <- c(numlist[2])
c <- c(colorlist[2])
pdf(file = "Res/20170721_numbersFinalPlotData/circo_all_psite_thisStudy.pdf", width = 7, height = 7)
# circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.barplot(l = l, col = c, bg.border = "white", track.height = 0.6)
dev.off()


# only protein this study
l <- c(numlist[1])
c <- c(colorlist[1])
pdf(file = "Res/20170721_numbersFinalPlotData/circo_all_protein_thisStudy.pdf", width = 7, height = 7)
# circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.barplot(l = l, col = c, bg.border = "white", track.height = 0.6)
dev.off()

# only this study, include phospho-site and proteins
l <- c(numlist)
c <- c(colorlist)
pdf(file = "Res/20170721_numbersFinalPlotData/circo_all_prot&site_thisStudy.pdf", width = 7, height = 7)
# circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.barplot(l = l, col = c, bg.border = "white", track.height = 0.6)
dev.off()

##################################### final figures #########################################
source("R/final/circos/circos.trackBarplot.R")

# only phosphosite this study
l <- c(numlist[2])
c <- c(colorlist[2])
tiff(file = "pubFigures/figure1_psite_numbers.tiff", width = 6, height = 6, units = "cm", res = 300)
circos.barplot(l = l, col = c, bg.border = "white", track.height = 0.6)
dev.off()


# only protein this study
l <- c(numlist[1])
c <- c(colorlist[1])
tiff(file = "pubFigures/figure1_prot_numbers.tiff", width = 6, height = 6, units = "cm", res = 300)
circos.barplot(l = l, col = c, bg.border = "white", track.height = 0.6, axis.sep = 2000, 
               axis.step = 1000, pat = "000", pat.rep = "K")
dev.off()
