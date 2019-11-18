library(scatterD3)
library(metap)

# function def
proclist <- function(x) {
  x <- x[!is.na(x$nobs), ]
  x <- x[x$nobs > 3, ]
  p1side <- x$pval/2
  p1side[x$r < 0] <- 1 - p1side[x$r < 0]
  x$pval1sided <- p1side
  x
}

linearScale <- function(x, min, max) {
  x <- x - min(x)
  x * (max-min)/max(x)+ min
}
transcolor <- function(x, alpha = 1) {
  if (alpha <= 1)
    alpha <- 255*alpha
  x <- col2rgb(x)
  rgb(x[1, ], x[2,], x[3, ], alpha = alpha, maxColorValue = 255)
}

complexp <- function(a, oneSide = FALSE) {
  xa <- a$pval
  if (oneSide) 
    xa <- a$pval1sided
  r <- tapply(xa, a$complex, function(x) {
    # r <- min(x)
    if (length(x) == 1) {
      r <- x
    } else if (all(x == 1)) {
      r <- 1
    } else {
      # x[x == 1] <- 1-1e-16
      # x[x == 0] <- 1e-16
      r <- wilkinsonp(x)$p
    }
    # -log10(r)
    c(min(16, -log10(r)), min(16, -log10(min(x))))
  })
  df <- as.data.frame(do.call(rbind, r))
  colnames(df) <- c("combp", "minp")
  df
} 

savefigures <- function(a, b, fname) {
  d <- b
  pa <- complexp(a, oneSide = FALSE)
  pd <- complexp(d, oneSide = TRUE)
  table(pa[, 1] > 6)
  table(pd[, 1] > 6)
  it <- intersect(rownames(pa), rownames(pd))
  
  fn <- file.path("Res/20170927_corumCorrelationPlots", paste(fname, "cor.png", sep = "_"))
  png(fn, width = 20, height = 20, units = "cm", res = 300)
  size <- pmax(pa[it, 1], pa[it, 2])
  color <- rep("red", length(it))
  color[pa[it, 1] < pa[it, 2]] <- "green"
  plot(pa[it, 1], pd[it, 1], pch = 19, cex = linearScale(size, min = 1, max = 3), col = transcolor(color, 0.4), 
       frame.plot = FALSE, xlab = "-log10(p value) Phospho data ", ylab = "-log10(p value) Protein data ")
  axis(side = 1)
  axis(side = 2)
  abline(a = 0, b = 1, lty = 2)
  abline(h = -log10(0.05), v = -log10(0.05), lty = 3)
  dev.off()
  
  #
  fn <- file.path("Res/20170927_corumCorrelationPlots", paste(fname, "bar.png", sep = "_"))
  png(fn, width = 15, height = 20, units = "cm", res = 300)
  pu <- setdiff(rownames(pa), rownames(pd))
  lp <- structure(pa[pu, 1], names = pu)
  par(mar = c(20, 4, 1, 1))
  barplot(sort(lp, decreasing = TRUE)[1:min(15, length(lp))], las = 2, ylab = "-log10(p value) Phospho data")
  abline(h = -log10(0.05), lty = 2)
  dev.off()
  
  scatterD3(pa[it, 1], pd[it, 1], lab = it)
}

## loading data
l <- readRDS(file = "Res/20170926_corumCorrelation/extractResults.RDS")

a <- proclist(l$site.trypsin.nci60)
b <- proclist(l$site.gluc.nci60)
c <- proclist(l$site.trypsin.crc65)

d <- proclist(l$prot.trypsin.nci60)
e <- proclist(l$prot.gluc.nci60)
f <- proclist(l$prot.trypsin.crc65)

sl <- list(site.trypsin.nci60 = a, 
           site.gluc.nci60 = b, 
           site.trypsin.crc65 = c,
           prot.trypsin.nci60 = d,
           prot.gluc.nci60 = e,
           prot.trypsin.crc65 = f)

saveRDS(sl, file = "Res/20170926_corumCorrelation/summarizedRes.RDS")

savefigures(a = a, b = d, fname = "nci60.trypsin")
savefigures(a = b, b = e, fname = "nci60.gluc")
savefigures(a = c, b = f, fname = "crc65.trypsin")


