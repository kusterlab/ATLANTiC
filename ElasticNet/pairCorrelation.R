library(fastmatch)
library(scatterD3)

res <- readRDS(file = "Res/20170926_corumCorrelation/summarizedRes.RDS")

reducePairs <- function(a) {
  a0 <- paste(a$subunit1, a$subunit2, a$pval)
  ar <- apply(a, 2, function(x) {
    tapply(x, a0, function(y) {
      paste(unique(y), collapse = ";")
    })
  })
  ar <- as.data.frame(ar, stringsAsFactors = FALSE)
  ar$pval <-as.numeric(ar$pval)
  ar$r <-as.numeric(ar$r)
  ar$rImpute <-as.numeric(ar$rImpute)
  ar$nobs <-as.numeric(ar$nobs)
  ar$pval1sided <-as.numeric(ar$pval1sided)
  ar$pair <- paste(ar$subunit1, ar$subunit2)
  ar
}

phosvsprot <- function(a, b) {
  # a phos
  # b prot
  
  a0 <- reducePairs(a)
  b0 <- reducePairs(b)
  
  rownames(a0) <- a0$pair
  rownames(b0) <- b0$pair
  
  apair <- unique(c(a0$pair, b0$pair))
  a0 <- a0[apair, ]
  b0 <- b0[apair, ]
  
  colnames(a0) <- paste0("phos.", colnames(a0))
  colnames(b0) <- paste0("prot.", colnames(b0))
  
  ab <- cbind(a0, b0)
  rownames(ab) <- apair
  ab
}

#
nci60.tryp <- phosvsprot(res$site.trypsin.nci60, res$prot.trypsin.nci60)
nci60.gluc <- phosvsprot(res$site.gluc.nci60, res$prot.gluc.nci60)
crc65.tryp <- phosvsprot(res$site.trypsin.crc65, res$prot.trypsin.crc65)

#### ======================= CRC65 =====================
a <- pmin(-log10(crc65.tryp$phos.pval), 16)
a[is.na(a)] <- 0
b <- pmin(-log10(crc65.tryp$prot.pval), 16)
b[is.na(b)] <- 0


sum(a <3 & b > 6)
i <- a >6 & b < 3
f <- crc65.tryp[i, ]
scatterD3(a[i], b[i], lab = paste(f$phos.complex, f$phos.feature1, f$phos.feature2, sep = "_"))



png("Res/20170927_corumCorrelationPlots/complexPairCor_prot_vs_phos_CRC65.png", width = 10, height = 10, units = "cm", res = 300)
plot(a, b, xlab = "-log10(pvalue) site pairs", ylab = "-log10(pvalue) protein pairs", 
     pch = 20, cex = 0.6, col = "gray", main = "CRC65")
segments(6, 3, 6, -1, lty = 2)
segments(6, 3, 17, 3, lty = 2)
segments(3, 6, 3, 17, lty = 2)
segments(3, 6, -1, 6, lty = 2)
dev.off()




## ======================= NCI60 trypsin =====================
a <- pmin(-log10(nci60.tryp$phos.pval), 16)
a[is.na(a)] <- 0
b <- pmin(-log10(nci60.tryp$prot.pval), 16)
b[is.na(b)] <- 0

sum(a <3 & b > 6)
i <- a >6 & b < 3
f <- nci60.tryp[i, ]
scatterD3(a[i], b[i], lab = paste(f$phos.complex, f$phos.feature1, f$phos.feature2, sep = "_"))

png("Res/20170927_corumCorrelationPlots/complexPairCor_prot_vs_phos_NCI60.trypsin.png", width = 10, height = 10, units = "cm", res = 300)
plot(a, b, xlab = "-log10(pvalue) site pairs", ylab = "-log10(pvalue) protein pairs", 
     pch = 20, cex = 0.6, col = "gray", main = "NCI60 trypsin")
segments(6, 3, 6, -1, lty = 2)
segments(6, 3, 17, 3, lty = 2)
segments(3, 6, 3, 17, lty = 2)
segments(3, 6, -1, 6, lty = 2)
dev.off()





#### ======================= NCI60 gluc =====================
a <- pmin(-log10(nci60.gluc$phos.pval), 16)
a[is.na(a)] <- 0
b <- pmin(-log10(nci60.gluc$prot.pval), 16)
b[is.na(b)] <- 0

sum(a <3 & b > 6)
i <- a >6 & b < 3
f <- nci60.gluc[i, ]
scatterD3(a[i], b[i], lab = paste(f$phos.complex, f$phos.feature1, f$phos.feature2, sep = "_"))

png("Res/20170927_corumCorrelationPlots/complexPairCor_prot_vs_phos_NCI60.gluc.png", width = 10, height = 10, units = "cm", res = 300)
plot(a, b, xlab = "-log10(pvalue) site pairs", ylab = "-log10(pvalue) protein pairs", 
     pch = 20, cex = 0.6, col = "gray", main = "NCI60 gluc")
segments(6, 3, 6, -1, lty = 2)
segments(6, 3, 17, 3, lty = 2)
segments(3, 6, 3, 17, lty = 2)
segments(3, 6, -1, 6, lty = 2)
dev.off()
