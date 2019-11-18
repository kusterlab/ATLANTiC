library(ade4)
library(svd)
#

source("R/concordAnalysis/concord.R")
source("R/concordAnalysis/concord.sparse.raw.R")


# source("R/drugResponseProcess/extract.drugmat.nci60.R")
# fp.nci60 <- readRDS("Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
# try.nci60 <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
# glu.nci60 <- readRDS("Res/20170110_mq15processLOD2/sites.gluC.RDS")
# 
# dmat <- drug.nci60$dtp.gi50
# dmat[is.na(dmat)] <- 0
# 
# prot <- fp.nci60$fun_prepElnet(fp.nci60, celllines = "all")
# prot <- prot$imputed
# tryp <- try.nci60$fun_prepElnet(try.nci60, celllines = "all", panel = "NCI60", regsitesOnly = TRUE)
# tryp <- tryp$imputed
# gluc <- glu.nci60$fun_prepElnet(glu.nci60, celllines = 'all', regsitesOnly = TRUE)
# gluc <- gluc$imputed
# 
# Xmats <- list(prot = prot, tryp = tryp, gluc = gluc)
# Ymat <- dmat
# 
# save(list = c("Xmats", "Ymat"), file = "/media/msdata5/users_files/Chen/Projects/Concordance/R/Functions/example.data.RDA")


load("/media/msdata5/users_files/Chen/Projects/Concordance/R/Functions/example.data.RDA")

X = Xmats
Y = Ymat

source("R/concordAnalysis/concord.R")

##
#
r1 <- concord(X = Xmats, Y = Ymat, nf = 5, dmod = 1)
r2 <- concord(X = Xmats, Y = Ymat, nf = 5, dmod = 2)
r3 <- concord(X = Xmats, Y = Ymat, nf = 5, dmod = 3)
r4 <- concord(X = Xmats, Y = Ymat, nf = 5, dmod = 4)
# r5 <- concord(X = Xmats, Y = Ymat, nf = 5, dmod = 5)

#
cbind(checkCC(r1), 
      checkCC(r2), 
      checkCC(r3), 
      checkCC(r4))

all.equal(r1, r2)
all.equal(r1, r3)
all.equal(r1, r4)
# all.equal(r1, r5)

all.equal(r2, r3)
all.equal(r2, r4)
all.equal(r2, r5)

all.equal(r3, r4)
all.equal(r3, r5)

all.equal(r4, r5)

cor(r1$loading, r4$loading)



round(cor(r1$ys, r1$gls), 2)
round(cor(r2$ys, r2$gls), 2)
round(cor(r3$ys, r3$gls), 2)

checkCC(r3)
sum(r3$var)
sum(r3$S^2)
(sum(r2$var) + sum(r3$var))/2



round(crossprod(r1$yloading), 2)
sum(r5$var)
sum(r5$S^2)

diag(cor(r4$ys, r4$gls)) - diag(cor(r5$ys, r5$gls))


round(cor(r4$ys, r4$gls), 2)
round(cor(r5$ys, r5$gls), 2)







library(mogsa)
rr <- sconcord(Xmats, Ymat, k = c(0.1, 0.2, 0.2), nf = 6, deflat = "globalLoading")


cor(rr$scoreXall, rr$scoreY)


plot(rr$scoreXall)


colSums(rr$loadingXs != 0)

sum(sapply(Xmats, nrow) * c(0.1, 0.2, 0.2))






(rr$scoreXs)

barplot(diag(cor(rr$scoreY, rr$scoreXall)))


cor(r1$yloading[, 1:6], rr$loadingY)
cor(r2$yloading[, 1:6], rr$loadingY)
cor(r3$yloading[, 1:6], rr$loadingY)
cor(r4$yloading[, 1:6], rr$loadingY)
cor(r5$yloading[, 1:6], rr$loadingY)






crossprod(r1$loading[, 1:6], rr$loadingXs)
crossprod(r2$loading[, 1:6], rr$loadingXs)
crossprod(r3$loading[, 1:6], rr$loadingXs)
crossprod(r4$loading[, 1:6], rr$loadingXs)
crossprod(r5$loading[, 1:6], rr$loadingXs)


crossprod(r1$yloading[, 1:6], rr$loadingY)
crossprod(r2$yloading[, 1:6], rr$loadingY)
crossprod(r3$yloading[, 1:6], rr$loadingY)
crossprod(r4$yloading[, 1:6], rr$loadingY)
crossprod(r5$yloading[, 1:6], rr$loadingY)




cor(r1$ys[, 1], rr$scoreY[, 1])
cor(r1$ys[, 2], rr$scoreY[, 2])
cor(r1$ys[, 1], rr$scoreY[, 2])

head(r1$ys[, 1:6])


r1$ys[, 1:6]
dim(rr$scoreY)
dim(r1$ys)

crossprod(rr$loadingXs)




barplot(diag(cor(r2$gls, r2$ys)))
barplot(diag(cor(r1$gls, r1$ys)))
barplot(diag(cor(r3$gls, r3$ys)))
barplot(diag(cor(r4$gls, r4$ys)))
barplot(diag(cor(r5$gls, r5$ys)))

#
sum(crossprod(r4$gls[, 1:59], r4$ys[, 1:59])^2)
sum(r4$var)

sum(diag(crossprod(r2$gls[, 1:59], r2$ys[, 1:59]))^2)
sum(r2$var)

sum(diag(crossprod(r3$gls[, 1:59], r3$ys[, 1:59]))^2)
sum(r3$var)


r3$gls - r2$gls
r2$gls

cor(r2$loading)[1:10, 1:10]
cor(r3$loading)[1:10, 1:10]

plot(r2$loading[, 2], r3$loading[, 2])
r2$loading[, 2] - r3$loading[, 2]

r2$loading[, 1:59] - r3$loading[, 1:59]
r2$gls[, 1:59] - r3$gls[, 1:59]

plot(r2$gls[, 6], r3$gls[, 6])
cor(r2$gls[, 1:10], r3$gls[, 1:10])

# ######################################
# r <- r2
# library(scatterD3)
# scatterD3(r$yloading[, 1], r$yloading[, 2], lab = rownames(Ymat), col_var = drug.nci60$dtp.annot$MOA)
# scatterD3(r$loading[, 1], r$loading[, 2], col_var = r$i.feature, lab = unlist(lapply(Xmats, rownames)))
# 
# scatterD3(r$yloading[, 3], r$yloading[, 4], lab = rownames(Ymat), col_var = drug.nci60$dtp.annot$MOA)
# scatterD3(r$loading[, 3], r$loading[, 4], col_var = r$i.feature, lab = unlist(lapply(Xmats, rownames)))
# 
# scatterD3(r$yloading[, 5], r$yloading[, 6], lab = rownames(Ymat), col_var = drug.nci60$dtp.annot$MOA)
# scatterD3(r$loading[, 5], r$loading[, 6], col_var = r$i.feature, lab = unlist(lapply(Xmats, rownames)))

