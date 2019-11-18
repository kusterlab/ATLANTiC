library(randomForest)
library(rpart)
library(matrixStats)
library(nnet)
library(e1071)
library(beeswarm)
library(parallel)

dat <- readRDS(file = "Res/20170317_Fluorouracil/prep.dat.RDS")

df <- dat$lod2(dat$crc65.prot)
df$ABCC4.1 <- NULL

crc65.resp <- readRDS("fromMartin/5FU.rda")
nn <- paste(rownames(crc65.resp), "_CRC65", sep = "")
nn[nn == "HCT 116_CRC65"] <- "HCT116_CRC65"
nn[nn == "HT-29_CRC65"] <- "HT29_CRC65"
crc65.resp <- structure(crc65.resp[, 1], names = nn)

##
fu <- crc65.resp[!is.na(crc65.resp)]
barplot(sort(fu), las = 2)
length(fu)
dim(df)

##@@@@
fu2 <- (dat$fu19893 + dat$fu757036)/2
res <- readRDS("Res/20170317_Fluorouracil/RFModelsPredCorCRC65.RDS")
boxplot(res)
df.nci <- dat$lod2(dat$nci60.prot)


oobr <- c(sapply(3:7, function(i) sort(res[[i]], decreasing = TRUE)[1:5]))
mks <- c(sapply(3:7, function(i) names(sort(res[[i]], decreasing = TRUE)[1:5])))

am <- colnames(df)


map <- sapply(mks, function(x) {
  mk <- strsplit(x, "\\+")[[1]]
  r1 <- am %in% mk
  names(r1) <- am
  
  rf <- randomForest(fu ~ ., df[names(fu), mk])
  pred <- predict(rf, df.nci[, mk])
  
  identical(names(pred), names(fu2))
  
  r2 <- sapply(unique(dat$nci.xref$too.short), function(too) {
    cor(pred[dat$nci.xref$too.short == too], fu2[dat$nci.xref$too.short == too], use = "pair")
  })
  
  r3 <- cor(pred, fu2, use= "pair")
  c(r1, r2, all.too=r3)
})

map <- rbind(map, Roob = oobr)
write.table(round(map, digits = 3), 
            file = "Res/20170317_Fluorouracil/fluTab.txt", col.names = TRUE, 
            row.names = TRUE, sep = "\t", quote = FALSE)

## Some examples
library(beeswarm)

bdf <- as.data.frame(t(map[16:25, ]))
cc <- unique(dat$nci.xref$color)
names(cc) <- unique(dat$nci.xref$too.short)
cc <- cc[names(bdf)]
cc[length(cc)] <- "pink"

png("Res/20170317_Fluorouracil/cor.coef.fu.too.tiff", width = 16, height = 12, units = "cm", res = 200)
beeswarm(bdf, col = cc, pch = 19, corral = "gutter", ylab = "correlation coefficent (R)")
bxplot(bdf, add = TRUE)
dev.off()



sig <- "TYMS+RRM1+UCK1"
mk <- strsplit(sig, "\\+")[[1]]
rf <- randomForest(fu ~ ., df[names(fu), mk])
pred <- predict(rf, df.nci[, mk])
identical(names(pred), names(fu2))

png("Res/20170317_Fluorouracil/fu.co.predvsmeasured.3.tiff", width = 10, height = 10, units = "cm", res = 200)
plot(pred[dat$nci.xref$too.short == "CO"], fu2[dat$nci.xref$too.short == "CO"], 
     col = "orange3", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured -log10(GI50)", main = sig)
dev.off()

###
sig <- "TYMP+RRM1+UPP1+UCK2+UCK1"
mk <- strsplit(sig, "\\+")[[1]]
rf <- randomForest(fu ~ ., df[names(fu), mk])
pred <- predict(rf, df.nci[, mk])
identical(names(pred), names(fu2))

png("Res/20170317_Fluorouracil/fu.co.predvsmeasured.5.tiff", width = 10, height = 10, units = "cm", res = 200)
plot(pred[dat$nci.xref$too.short == "CO"], fu2[dat$nci.xref$too.short == "CO"], 
     col = "orange3", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured -log10(GI50)", main = sig)
dev.off()





