source("https://raw.githubusercontent.com/mengchen18/RFunctionCollection/master/outliers.R")
library(matrixStats)

croot <- "/media/kusterlab/internal_projects/active/NCI60_Phospho/martin"
ev_phosprots_nci60 <- readRDS(file.path(croot, 'dat/preproc/ev_phosprots_nci60.rds'))
ev_phosprots_crc65 <- readRDS(file.path(croot, 'dat/preproc/ev_phosprots_crc65.rds'))


matCRC <- ev_phosprots_crc65
matCRC <- matCRC[, setdiff(colnames(matCRC), "CoCM-1_CRC65")]
dim(matCRC)
matNCI <- ev_phosprots_nci60

##
tiff("Res/20180328_outlierAnalysis/filterThresh_crc65_Pprot.tiff", width = 18, height = 6, units = "cm", res = 300)
i <- findOutlier(matCRC, foldthresh = 5, pvalue = 0.1, reachLowBound = FALSE, window = 0.2)
dev.off()

x <- i$outlierMatrix
x <- t(apply(x, 1, function(x1) {
  x1[is.na(x1)] <- min(x1, na.rm = TRUE)-log10(2)
  x1
}))
rownames(x) <- NULL

tiff("Res/20180328_outlierAnalysis/heatmap_outlier_crc65_Pprotein.tiff", width = 16, height = 16, units = "cm", res = 150)
pheatmap(x, color = c("white", colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)),  fontsize_col = 7, scale = "row", 
         breaks = c(-8, seq(2, 8, length.out = 100)))
dev.off()

lss1 <- apply(i$outlierMatrix, 1, which.max)
lss1 <- data.frame(gene = rownames(i$outlierMatrix), 
                   cellline = colnames(i$outlierMatrix)[lss1])
lss1 <- lss1[order(lss1$gene), ]
lss1 <- lss1[order(lss1$cellline), ]

####


tiff("Res/20180328_outlierAnalysis/filterThresh_nci60_Pprot.tiff", width = 18, height = 6, units = "cm", res = 300)
i <- findOutlier(matNCI, foldthresh = 5, pvalue = 0.1, reachLowBound = FALSE, window = 0.2)
dev.off()

x <- i$outlierMatrix
x <- t(apply(x, 1, function(x1) {
  x1[is.na(x1)] <- min(x1, na.rm = TRUE)-log10(2)
  x1
}))
rownames(x) <- NULL

tiff("Res/20180328_outlierAnalysis/heatmap_outlier_nci60_Pprotein.tiff", width = 16, height = 16, units = "cm", res = 150)
pheatmap(x, color = c("white", colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)),  fontsize_col = 7, scale = "row", 
         breaks = c(-8, seq(2, 8, length.out = 100)))
dev.off()

lss2 <- apply(i$outlierMatrix, 1, which.max)
lss2 <- data.frame(gene = rownames(i$outlierMatrix), 
                   cellline = colnames(i$outlierMatrix)[lss2])
lss2 <- lss2[order(lss2$gene), ]
lss2 <- lss2[order(lss2$cellline), ]

lss1$panel <- "CRC65"
lss2$panel <- "NCI60"
lss <- rbind(lss1, lss2)
write.table(lss, "Res/20180328_outlierAnalysis/outlierPprot.txt", col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = "\t")
