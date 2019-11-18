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

df <- dat$lod2(dat$nci60.prot)
df <- cbind(df, dat$lod2(dat$nci60.tk1phos))
df$TK1 <- NULL
df$UMPS.1 <- NULL

# lm model
res <- lapply(1:15, function(n) {
  print(n)
  combs <- combn(colnames(df), n, simplify = FALSE)
  sigs <- sapply(combs, paste, collapse = "+")

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

saveRDS(res, file = "Res/20170317_Fluorouracil/linearModelsPredCorTK1phos.RDS")

## Random forest
res <- lapply(1:15, function(n) {
  print(n)
  combs <- combn(colnames(df), n, simplify = FALSE)
  sigs <- sapply(combs, paste, collapse = "+")
  
  r <- mclapply(combs, function(g) {
    dff <- df[, g, drop = FALSE]
    i <- !is.na(fu)
    rf <- randomForest(fu[i] ~ ., dff[i, , drop = FALSE])
    min(rf$mse)
  }, mc.cores = 20)
  r <- unlist(r)
  names(r) <- sigs
  r
})

saveRDS(res, "Res/20170317_Fluorouracil/RFModelsPredCorTK1phos.RDS")
