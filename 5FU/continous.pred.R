library(randomForest)
library(rpart)
library(matrixStats)
library(nnet)
library(e1071)
library(beeswarm)
library(parallel)

dat <- readRDS(file = "Res/20170317_Fluorouracil/prep.dat.RDS")
plot(dat$fu19893, dat$fu757036)

### Continous predictive model
# use mean
fu <- (dat$fu19893 + dat$fu757036)/2

# lm model
res <- lapply(1:15, function(n) {
  print(n)
  combs <- combn(colnames(dat$nci60.prot), n, simplify = FALSE)
  sigs <- sapply(combs, paste, collapse = "+")
  df <- dat$lod2(dat$nci60.prot)
  
  r <- mclapply(combs, function(g) {
    dff <- df[, g, drop = FALSE]
    cor <- replicate(100, {
      ts <- sample(names(fu), replace = TRUE)
      es <- setdiff(names(fu), ts)
      lm <- lm(fu[ts] ~ ., dff[ts, , drop = FALSE])
      v <- predict(lm, dff[es, , drop = FALSE])
      cor(fu[es], v, use = "pair")
    }, simplify = TRUE)
  }, mc.cores = 20)
  r <- do.call("cbind", r)
  colnames(r) <- sigs
  r
})

saveRDS(res, file = "Res/20170317_Fluorouracil/linearModelsPredCor.RDS")
res2 <- readRDS(file = "Res/20170317_Fluorouracil/linearModelsPredCor.RDS")
rr <- do.call(cbind, res2)
cmr <- colMedians(rr)
names(cmr) <- colnames(rr)

hist(cmr, breaks = 100)
sort(cmr, decreasing = TRUE)[1:10]

plot(fu, dat$lod2(dat$nci60.prot)[, "UCK2"])
cor(fu, dat$lod2(dat$nci60.prot)[, "UCK2"], use = "pair")

## Random forest

res <- lapply(1:15, function(n) {
  print(n)
  combs <- combn(colnames(dat$nci60.prot), n, simplify = FALSE)
  sigs <- sapply(combs, paste, collapse = "+")
  df <- dat$lod2(dat$nci60.prot)
  
  r <- mclapply(combs, function(g) {
    dff <- df[, g, drop = FALSE]
    i <- !is.na(fu)
    rf <- randomForest(fu[i] ~ ., dff[i, ])
    min(rf$mse)
  }, mc.cores = 20)
  r <- unlist(r)
  names(r) <- sigs
  r
})

saveRDS(res, "Res/20170317_Fluorouracil/RFModelsPredCor.RDS")

#
boxplot(res[-1])
which.min(res[[3]])

i <- !is.na(fu)
df <- dat$lod2(dat$nci60.prot)

crc65.resp <- readRDS("fromMartin/5FU.rda")
nn <- paste(rownames(crc65.resp), "_CRC65", sep = "")
nn[nn == "HCT 116_CRC65"] <- "HCT116_CRC65"
nn[nn == "HT-29_CRC65"] <- "HT29_CRC65"
crc65.resp <- structure(crc65.resp[, 1], names = nn)

mks <- names(which.min(res[[7]]))
mks <- strsplit(mks, "\\+")[[1]]
df.crc <- dat$lod2(dat$crc65.prot)

mks <- c("PPAT", "UCK2")
###################################
rf <- randomForest(fu[i] ~ ., df[i, mks])
pred <- predict(rf, df.crc[, mks])
###################################
lmod <- lm(fu[i] ~ ., df[i, mks])
pred <- predict(lmod, df.crc[, mks])

###
plot(pred, crc65.resp[names(pred)])
cor.test(pred, crc65.resp[names(pred)], use = "pair")









