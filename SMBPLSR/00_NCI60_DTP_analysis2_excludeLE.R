library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)
library(omic3plus)
library(impute)
library(WebGestaltR)

dat <- readRDS("Res/20170904_concordanceData/nci60_dtp.RDS")

icell <- dat$xref$too.short != "LE"
y <- dat$y$data[, icell]
X <- lapply(dat$x, function(x) x$imputed[, icell])

y <- y[rowMaxs(y) > 1.5, ]
X <- lapply(X, function(x) {
  x[rowSds(x)>1e-5, ]
})

sapply(X, dim)


cho.ky <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
cho.kx <- 50*(2^(0:3))

vp <- concord(X, y, kx = cho.kx, ky = cho.ky, option = "nk", ncomp = 40,
              scan = FALSE, ncores = 5,
              dmod = 1, center.x = TRUE, scale.x = TRUE,
              center.y = TRUE, scale.y = TRUE, pos = FALSE)

loadingy <- data.frame(vp$loading.y)
loadingy <- cbind(drug = rownames(y), dat$drugannot[rownames(y), "MOA"], loadingy)

####

st <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
sg <- readRDS("Res/20170110_mq15processLOD2/sites.gluC.RDS")
pt <- readRDS("Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
pg <- readRDS("Res/20170719_nci60.proteinGroups.gluc/nci60.proteinGroups.gluc.RDS")


wb <- createWorkbook("res")
addWorksheet(wb, "DrugLoadings")
writeDataTable(wb = wb, x = loadingy, sheet = 1)

for (i in 1:ncol(vp$score.x)) {
  print(i)
  ir <- abs(vp$loading.x[, i]) > 1e-7
  a <- vp$loading.x[ir, i]
  a <- data.frame(id = names(a), loading = a, dataset = vp$loading.x.index[ir])
  a <- a[order(a$loading), ]
  a$symbol <- sapply(strsplit(rownames(a), split = " ", fixed = 2), "[", 2)
  
  n <- rownames(a)
  ast <- st$annot[intersect(n, rownames(st$annot)), c("Protein names", "Position", "Known site", "Sequence window", "Regulatory site")]
  asg <- sg$annot[intersect(n, rownames(sg$annot)), c("Protein names", "Position", "Known site", "Sequence window", "Regulatory site")]
  apt <- pt$annot[intersect(n, rownames(pt$annot)), c("Protein names"), drop = FALSE]
  apg <- pg$annot[intersect(n, rownames(pg$annot)), c("Protein names"), drop = FALSE]
  if (nrow(apt) > 0)
    apt$"Position" <- apt$"Known site" <- apt$"Sequence window" <- apt$"Regulatory site" <- NA
  if (nrow(apg) > 0)
    apg$"Position" <- apg$"Known site" <- apg$"Sequence window" <- apg$"Regulatory site" <- NA
  ar <- rbind(ast, asg, apt, apg)
  m <- cbind(a, ar[n, ])
  
  addWorksheet(wb, sheetName = paste0("comp", i))
  writeDataTable(wb = wb, x = m, sheet = 1+i)
}

saveWorkbook(wb, "Res/20170915_concordResSummary/concord-c40.xlsx", overwrite = TRUE)

i <- 9
plot(vp$score.x[, i], vp$score.y[, i], col = dat$xref$color[icell], pch = 20)
saveRDS(vp, "Res/20170915_concordResSummary/concord_noLE_comp40.RDS")

