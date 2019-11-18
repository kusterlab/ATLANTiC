library(randomForest)
library(rpart)
library(matrixStats)
library(nnet)
library(e1071)
library(beeswarm)
library(parallel)

dat <- readRDS(file = "Res/20170317_Fluorouracil/prep.dat.RDS")
df <- dat$lod2(dat$crc65.prot)
df <- cbind(df, dat$lod2(dat$crc65.tk1phos))
df$ABCC4.1 <- NULL
df$TK1 <- NULL
df$UMPS.1 <- NULL

crc65.resp <- readRDS("fromMartin/5FU.rda")
nn <- paste(rownames(crc65.resp), "_CRC65", sep = "")
nn[nn == "HCT 116_CRC65"] <- "HCT116_CRC65"
nn[nn == "HT-29_CRC65"] <- "HT29_CRC65"
crc65.resp <- structure(crc65.resp[, 1], names = nn)

##
fu <- crc65.resp[!is.na(crc65.resp)]
barplot(sort(fu))
length(fu)
dim(df)


# # lm model
# res <- lapply(1:15, function(n) {
#   print(n)
#   combs <- combn(colnames(df), n, simplify = FALSE)
#   sigs <- sapply(combs, paste, collapse = "+")
#   
#   r <- mclapply(combs, function(g) {
#     dff <- df[, g, drop = FALSE]
#     cor <- replicate(100, {
#       ts <- sample(names(fu), replace = TRUE)
#       es <- setdiff(names(fu), ts)
#       lm <- lm(fu[ts] ~ ., dff[ts, , drop = FALSE])
#       v <- predict(lm, dff[es, , drop = FALSE])
#       cor(fu[es], v, use = "pair")
#     }, simplify = TRUE)
#   }, mc.cores = 20)
#   r <- do.call("cbind", r)
#   colnames(r) <- sigs
#   r
# })

# saveRDS(res, file = "Res/20170317_Fluorouracil/linearModelsPredCorCRC65.RDS")

## Random forest
res <- lapply(1:15, function(n) {
  print(n)
  combs <- combn(colnames(df), n, simplify = FALSE)
  sigs <- sapply(combs, paste, collapse = "+")
  
  r <- mclapply(combs, function(g) {
    dff <- df[names(fu), g, drop = FALSE]
    rf <- randomForest(fu ~ ., dff, ntree = 10)
    sqrt(max(rf$rsq, na.rm = TRUE))
  }, mc.cores = 20)
  r <- unlist(r)
  names(r) <- sigs
  r
})

saveRDS(res, "Res/20170317_Fluorouracil/RFModelsPredCorCRC65TK1phos.RDS")

###
boxplot(res)

df.nci <- dat$lod2(dat$nci60.prot)

mks <- names(which.max(res[[3]]))
mks <- strsplit(mks, "\\+")[[1]]
rf <- randomForest(fu ~ ., df[names(fu), mks])
pred <- predict(rf, df.nci[, mks])
identical(names(pred), names(fu2))

plot(pred, fu2)
cor(pred, fu2, use = "pair")





