library(matrixStats)

tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
crc.phos <- readRDS("Res/20170721_numbersFinalPlotData/numberid_crc65_PHOS.RDS")
crc.fp <- readRDS("Res/20170721_numbersFinalPlotData/numberid_crc65_FP.RDS")
nci.phos <- readRDS("Res/20170721_numbersFinalPlotData/numberid_nci60_PHOS.RDS")
nci.fp <- readRDS("Res/20170721_numbersFinalPlotData/numberid_nci60_FP.RDS")


transparentColor <- function(x, alpha = 100) {
  rgbs <- col2rgb(x)
  rgb(rgbs[1, ], rgbs[2, ], rgbs[3, ], alpha = 100, maxColorValue = 255)
}

# NCI60 
# phos
xref <- tryp$nci60_xref
xref <- xref[order(xref$too.short), ]

png("Res/20170721_numbersFinalPlotData/numbers.nci60.phos.png", width = 18, height = 12, unit = "cm", res = 300)
cc <- transparentColor(unique(xref$color))
bp <- barplot(nci.phos$tissue["tisid", ], width = nci.phos$tissue["nsample", ], 
              col = cc, space = 0, border = NA, names.arg = "")
mtext(unique(xref$too.short), side = 3, line = 0, at = c(bp), col = unique(xref$color))
barplot(nci.phos$cellline, las = 2, add = TRUE, space = 0, axes = FALSE, legend.text = c("Both", "Trypsin", "GluC"),
        args.legend = list(x = 8, y = 11000, cex = 0.7), 
        cex.names = 0.7, names.arg = gsub("_NCI60", "", colnames(nci.phos$cellline)), las = 2)

med <- round(median(colSums(nci.phos$cellline)))
abline(h = med, lty = 3)
mtext(med, at = med+0.03*max(nci.phos$tissue), side = 4, las = 2, line = -0.6)
dev.off()


# prot
xref <- tryp$nci60_xref
xref <- xref[order(xref$too.short), ]

png("Res/20170721_numbersFinalPlotData/numbers.nci60.prot.png", width = 18, height = 12, unit = "cm", res = 300)
cc <- transparentColor(unique(xref$color))
bp <- barplot(nci.fp$tissue["tisid", ], width = nci.fp$tissue["nsample", ], ylim = c(0, 10000),
              col = cc, space = 0, border = NA, names.arg = "")
mtext(unique(xref$too.short), side = 3, line = 0, at = c(bp), col = unique(xref$color))
barplot(nci.fp$cellline, las = 2, add = TRUE, space = 0, axes = FALSE, legend.text = c("Both", "Trypsin", "GluC"),
        args.legend = list(x = 8, y = 3000, cex = 0.7), 
        cex.names = 0.7, names.arg = gsub("_NCI60", "", colnames(nci.fp$cellline)), las = 2)

med <- round(median(colSums(nci.fp$cellline)))
abline(h = med, lty = 3)
mtext(med, at = med+0.03*max(nci.fp$tissue), side = 4, las = 2, line = -0.6)
dev.off()

# CRC65 & NCI60
# protein
xref <- tryp$nci60_xref
xref <- xref[order(xref$too.short), ]
ccol <- unique(xref$color[xref$too.short == "CO"])

cellline <- cbind(nci.fp$cellline, rbind(0, crc.fp$cellline, 0))
wid <- c(nci.fp$tissue["nsample", ], crc.fp$tissue["nsample", ])
fcc <- c(unique(xref$color), ccol)
cc <- transparentColor(fcc)
tissue <- c(nci.fp$tissue["tisid", ], CRC65 = unname(crc.fp$tissue["tisid", ]))

png("Res/20170721_numbersFinalPlotData/numbers.nci60crc65.prot.png", width = 32, height = 12, unit = "cm", res = 300)
bp <- barplot(tissue, width = wid, ylim = c(0, 12000),
              col = cc, space = 0, border = NA, names.arg = "")
mtext(names(tissue), side = 3, line = 0, at = c(bp), col = fcc)
barplot(cellline, las = 2, add = TRUE, space = 0, axes = FALSE, legend.text = c("Both", "Trypsin", "GluC"),
        args.legend = list(x = 9, y = 3000, cex = 0.7), 
        cex.names = 0.6, names.arg = gsub("_NCI60|_CRC65|_GluC", "", colnames(cellline)), las = 2)
# med <- round(median(colSums(nci.fp$cellline)))
# abline(h = med, lty = 3)
# mtext(med, at = med+0.03*max(nci.fp$tissue), side = 4, las = 2, line = -0.6)
dev.off()

library(openxlsx)

wk <- createWorkbook()
addWorksheet(wk, "CELLLINE_NCI60_CRC65_protein")
addWorksheet(wk, "CELLLINE_NCI60_CRC65_phos")
addWorksheet(wk, "CELLLINE_NCI60_phos")
addWorksheet(wk, "CELLLINE_NCI60_protein")

writeData(wk, "CELLLINE_NCI60_CRC65_protein", x = t(cellline), colNames = TRUE, rowNames = TRUE)




# CRC65 & NCI60
# phos
xref <- tryp$nci60_xref
xref <- xref[order(xref$too.short), ]
ccol <- unique(xref$color[xref$too.short == "CO"])

cellline <- cbind(nci.phos$cellline, rbind(0, crc.phos$cellline, 0))
wid <- c(nci.phos$tissue["nsample", ], crc.phos$tissue["nsample", ])
fcc <- c(unique(xref$color), ccol)
cc <- transparentColor(fcc)
tissue <- c(nci.phos$tissue["tisid", ], CRC65 = unname(crc.phos$tissue["tisid", ]))

png("Res/20170721_numbersFinalPlotData/numbers.nci60crc65.phos.png", width = 32, height = 12, unit = "cm", res = 300)
bp <- barplot(tissue, width = wid, ylim = c(0, 50000),
              col = cc, space = 0, border = NA, names.arg = "")
mtext(names(tissue), side = 3, line = 0, at = c(bp), col = fcc)
barplot(cellline, las = 2, add = TRUE, space = 0, axes = FALSE, legend.text = c("Both", "Trypsin", "GluC"),
        args.legend = list(x = 9, y = 48000, cex = 0.7), 
        cex.names = 0.6, names.arg = gsub("_NCI60|_CRC65|_GluC", "", colnames(cellline)), las = 2)
# med <- round(median(colSums(nci.fp$cellline)))
# abline(h = med, lty = 3)
# mtext(med, at = med+0.03*max(nci.fp$tissue), side = 4, las = 2, line = -0.6)
dev.off()



writeData(wk, "CELLLINE_NCI60_CRC65_phos", x = t(cellline), colNames = TRUE, rowNames = TRUE)
writeData(wk, "CELLLINE_NCI60_phos", x = t(nci.phos$cellline), colNames = TRUE, rowNames = TRUE)
writeData(wk, "CELLLINE_NCI60_protein", x = t(nci.fp$cellline), colNames = TRUE, rowNames = TRUE)

saveWorkbook(wk, file = "Res/20170721_numbersFinalPlotData/numbers.xlsx")



