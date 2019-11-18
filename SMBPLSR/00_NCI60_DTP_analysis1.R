library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)
library(omic3plus)
library(impute)
library(WebGestaltR)
library(stringr)

dat <- readRDS("Res/20170904_concordanceData/nci60_dtp.RDS")
X <- lapply(dat$x, function(x) x$imputed)

cho.ky <- c(2, 4, 8, 16)
cho.kx <- c(50, 100, 200, 400, 800)

v <- concord(X, dat$y$data, kx = cho.kx, ky = cho.ky, option = "nk", fold = 5, ncores = 5, ncomp = 3, scan = FALSE,
             dmod = 1, center.x = TRUE, scale.x = FALSE, center.y = TRUE, scale.y = FALSE, pos = FALSE)
vp <- v

i <- 2
plot(vp$score.x[, i], vp$score.y[, i], col = dat$xref$color, pch = 20)
abline(a = 0, b = 1)
sel <- ext(vp, n = i, vy = rownames(dat$y$data))
selx <- sel$f
sely <- sel$d


gs <- str_split_fixed(selx$name, " |;", 3)
ref <- str_split_fixed(unlist(lapply(dat$x, function(x) rownames(x$imputed))), " |;", 3)


s <- WebGestaltR(enrichMethod = "ORA", 
                 organism = "hsapiens",
                 enrichDatabase="geneontology_Biological_Process",
                 minNum=3, maxNum=500,
                 interestGene = gs[, 2],
                 interestGeneType = "genesymbol",
                 referenceGene = ref[, 2], 
                 referenceGeneType = "genesymbol",
                 is.output = FALSE
)


##########################
ext <- function(obj, n = 1, vy) {
  d <- data.frame(name = vy[obj$loading.y[, n] != 0],
                  coef = obj$loading.y[obj$loading.y[, n] != 0, n], 
                  direction = 1)
  d$direction[d$coef < 0] <- -1
  d <- cbind(d, dat$drugannot[as.character(d$name), ])
  
  f <- data.frame(name = rownames(obj$loading.x)[obj$loading.x[, n] != 0], 
                  coef = obj$loading.x[obj$loading.x[, n] != 0, n], 
                  data = obj$loading.x.index[obj$loading.x[, n] != 0],
                  direction = 1)
  f$direction[f$coef < 0] <- -1
  list(f = f, d = d)
}


