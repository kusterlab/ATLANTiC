setwd("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen")
source("https://raw.githubusercontent.com/mengchen18/RFunctionCollection/master/outliers.R")

slist <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")

# library(stringr)
# library(pheatmap)
library(matrixStats)


if (!dir.exists("Res/20180328_outlierAnalysis"))
  dir.create("Res/20180328_outlierAnalysis")

matCRC65 <- slist$crc65_expr
matCRC65 <- matCRC65[, setdiff(colnames(matCRC65), "CoCM-1_CRC65")]
matNCI60 <- slist$nci60_expr


tiff("Res/20180328_outlierAnalysis/filterThresh_crc65.tiff", width = 18, height = 6, units = "cm", res = 300)
iCRC65 <- findOutlier(matCRC65, foldthresh = 5, pvalue = 0.1)
dev.off()

tiff("Res/20180328_outlierAnalysis/filterThresh_nci60.tiff", width = 18, height = 6, units = "cm", res = 300)
iNCI60 <- findOutlier(matNCI60, foldthresh = 5, pvalue = 0.1)
dev.off()

isx <- c(iCRC65$outlierIndexColumns, iNCI60$outlierIndexColumns)

annot <- slist$annot
colnames(annot)

retainCols <- c("Gene name", "Position", "Sequence window", "Regulatory site", "Pfam domains")


ls <- lapply(names(isx), function(x) {
  nm <- strsplit(x, "_")[[1]]
  print(nm)
  an <- annot[isx[[x]], retainCols]
  an$cell <- nm[1]
  an$panel <- nm[2]
  an
})


tab <- do.call(rbind, ls)
tab$AA <- substr(tab$`Sequence window`, 16, 16)


write.table(tab, file = "Res/20180328_outlierAnalysis/outlierTable.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
saveRDS(iCRC65, file = "Res/20180328_outlierAnalysis/outlier.CRC65.RDS")
saveRDS(iNCI60, file = "Res/20180328_outlierAnalysis/outlier.NCI60.RDS")



library(pheatmap)
library(RColorBrewer)
x <- iCRC65$outlierMatrix
x[is.na(x)] <- 0
rownames(x) <- NULL

tiff("Res/20180328_outlierAnalysis/heatmap_outlier_crc65.tiff", 
     width = 16, height = 16, units = "cm", res = 150)
pheatmap(x, color = c("white", colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)),  fontsize_col = 7)
dev.off()


x <- iNCI60$outlierMatrix
x[is.na(x)] <- 0
rownames(x) <- NULL
tiff("Res/20180328_outlierAnalysis/heatmap_outlier_nci60.tiff", 
     width = 16, height = 16, units = "cm", res = 150)
pheatmap(x, color = c("white", colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)),  fontsize_col = 7)
dev.off()




