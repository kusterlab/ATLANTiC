library(XML)
library(matrixStats)
dtp2 <- readHTMLTable("Dat/DTP_20170328/DTP_NCI60_ZSCORE.html", header = TRUE)
dtp2 <- dtp2[[1]]
dmat <- dtp2[, 7:66]
dmat <- apply(dmat, 2, function(x) as.numeric(as.character(x)))
dannot <- dtp2[, c(1:6, 67:68)]
rownames(dmat) <- rownames(dannot)

###
isens <- lapply(1:ncol(dmat), function(i) {
  which(dmat[, i] - rowMaxs(dmat[, -i], na.rm = TRUE) > 3)
})
names(isens) <-  colnames(dmat)
sort(sapply(isens, length))

qq <- lapply(names(isens), function(x) {
  xx <- isens[[x]]
  tm <- dannot[xx, ]
  tm$cellline <- x
  tm
})
qq <- do.call(rbind, qq)
qq <- cbind(qq, dmat[rownames(qq), ])

write.table(qq, file = "Res/20180328_outlierAnalysis/outlierDTP.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")




###################### beeswarms plot =================================
library(beeswarm)

nci <- readRDS("Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
xref <- nci$xref
xref <- xref[xref$name_DTP != "", ]

xmap <- read.delim("Res/20180328_outlierAnalysis/xref2.txt", header = FALSE, stringsAsFactors = FALSE)
identical(xmap$V1, xref$name_MQ)

##
qq0 <- qq[nchar(as.character(qq$`Mechanism of action`)) > 1, ]
qq0$`Drug name`

alki <- list(ldk378_776422 = dmat[which(dannot$`NSC #` == "776422"), xmap$V2], # LDK-378
             ldk378_771193 = dmat[which(dannot$`NSC #` == "777193"), xmap$V2], # LDK-378
             crizotinib = dmat[which(dannot$`Drug name` == "Crizotinib"), xmap$V2])

abli <- list(imatinib = dmat[which(dannot$`Drug name` == "Imatinib"), xmap$V2], 
             nilotinib = dmat[which(dannot$`Drug name` == "Nilotinib"), xmap$V2])

parpi <- list(olaparib_753686 = dmat[which(dannot$`NSC #` == "753686"), xmap$V2])

lnci60 <- c(alki, abli, parpi)
transcolor <- function(x) {
  cb <- col2rgb(x)
  rgb(cb[1, ], cb[2, ], cb[3, ], alpha  = 150, maxColorValue = 255)
}
cc <- replicate(6, transcolor(xref$color), simplify = FALSE)

tiff("Res/20180328_outlierAnalysis/nci60_drug_outliers_beeswarm.tiff", width = 12, height = 8, units = "cm", res = 150)
par(mar = c(6, 2, 0.5, 0.5))
boxplot(lnci60, border = "white", frame.plot = FALSE, pch = NA, las = 2)
beeswarm(lnci60, las = 2, corral = "wrap", las = 2, pwcol = cc, pch = 19, add = TRUE)
dev.off()





