source("https://raw.githubusercontent.com/mengchen18/RFunctionCollection/master/outliers.R")

library(matrixStats)
if (!dir.exists("Res/20180328_outlierAnalysis"))
  dir.create("Res/20180328_outlierAnalysis")


crc <- readRDS("Res/20170110_mq15processLOD2/crc65.proteinGroups.RDS")
matCRC <- crc$ibaq
matCRC <- matCRC[, setdiff(colnames(matCRC), "CoCM-1_CRC65")]
dim(matCRC)

nci <- readRDS("Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
matNCI <- nci$ibaq
dim(matNCI)


getab <- function(m, a, col = c("Gene names", "iBAQ", "Majority protein IDs", "Protein names")) {
  i <- findOutlier(m, foldthresh = 5, pvalue = 0.1, reachLowBound = FALSE)
  ls <- lapply(names(i$outlierIndexColumns), function(ii) {
    ir <- i$outlierIndexColumns[[ii]]
    d <- a[ir, col]
    if (nrow(d) == 0)
      return(d)
    d$cellline <- ii
    d
  })
  do.call(rbind, ls)
}

tiff("Res/20180328_outlierAnalysis/filterThresh_crc65_prot.tiff", width = 18, height = 6, units = "cm", res = 300)
ls1 <- getab(m = matCRC, a = crc$annot)
dev.off()

tiff("Res/20180328_outlierAnalysis/filterThresh_nci60_prot.tiff", width = 18, height = 6, units = "cm", res = 300)
ls2 <- getab(m = matNCI, a = nci$annot)
dev.off()


lss <- rbind(ls1, ls2)

write.table(lss, file = "Res/20180328_outlierAnalysis/outlierTableProtein.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

####
# mm <- nci$ibaq[which(nci$annot$`Gene name` == "ABCB1"), ]
mm <- nci$ibaq[which(nci$annot$`Gene name` == "PITX2"), ]
par(mar = c(10, 4, 1, 1))
barplot(sort(10^mm), las = 2)




####
i <- findOutlier(matCRC, foldthresh = 5, pvalue = 0.1, reachLowBound = FALSE)
x <- i$outlierMatrix
x <- t(apply(x, 1, function(x1) {
  x1[is.na(x1)] <- min(x1, na.rm = TRUE)-log10(2)
  x1
}))
# x[is.na(x)] <- 0
rownames(x) <- NULL

tiff("Res/20180328_outlierAnalysis/heatmap_outlier_crc65_protein.tiff", width = 16, height = 16, units = "cm", res = 150)
pheatmap(x, color = c("white", colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)),  fontsize_col = 7, scale = "row", 
         breaks = c(-8, seq(2, 8, length.out = 100)))
dev.off()

i <- findOutlier(matNCI, foldthresh = 5, pvalue = 0.1, reachLowBound = FALSE)
x <- i$outlierMatrix
x <- t(apply(x, 1, function(x1) {
  x1[is.na(x1)] <- min(x1, na.rm = TRUE)-log10(2)
  x1
}))
# x[is.na(x)] <- 0
rownames(x) <- NULL

tiff("Res/20180328_outlierAnalysis/heatmap_outlier_nci60_protein.tiff", width = 16, height = 16, units = "cm", res = 150)
pheatmap(x, color = c("white", colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)),  fontsize_col = 7, scale = "row", 
         breaks = c(-8, seq(2, 8, length.out = 100)))
dev.off()
