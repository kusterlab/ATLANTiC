library(randomForest)
library(rpart)
library(matrixStats)
library(nnet)
library(e1071)
library(beeswarm)
library(parallel)
library(rcellminer)
library(rcellminerData)
library(matrixStats)


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

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# fu2 <- (dat$fu19893 + dat$fu757036)/2
v <- getDrugActivityRepeatData("19893", onlyCellMinerExps = TRUE, concFormat = "NegLogGI50M")
dim(v)
cm <- cor(t(v), use = "pair")
hist(cm, breaks = 100)

boxplot(v, las = 2)
cm <- colMedians(v, na.rm = TRUE)
names(cm) <- colnames(v)
gi50um <- 10^(-cm) * 10^6
namevec <- structure(dat$nci.xref$name_MQ, names = dat$nci.xref$name_DTP)
gi50um <- gi50um[names(namevec)]
names(gi50um) <- namevec[names(gi50um)]
names(gi50um)[is.na(names(gi50um))] <- "MDAMB468_NCI60"

barplot(sort(gi50um[dat$nci.xref$too.short == "CO"]), las = 2)
# gi50um <- na.omit(gi50um)

barplot(sort(gi50um))
fu2 <- gi50um

## NCI data matrix
res <- readRDS("Res/20170317_Fluorouracil/RFModelsPredCorCRC65.RDS")
png("Res/20170324_Fluororacil_TCGA/markers.comb.boxplot.png", width = 18, height = 12, res = 200, units = "cm")
boxplot(res, xlab = "# of markers", ylab = "pseudo-R")
dev.off()

df.nci <- dat$lod2(dat$nci60.prot)

####### combinations and markers
oobr <- c(sapply(3:7, function(i) sort(res[[i]], decreasing = TRUE)[1:5]))
mks <- c(sapply(3:7, function(i) names(sort(res[[i]], decreasing = TRUE)[1:5])))

am <- colnames(df)
rf.all <- randomForest(fu ~ ., df[names(fu), ], importance = TRUE)
plot(rf.all)

png("Res/20170317_Fluorouracil/absScale/varImportance.absScale.png", 
    width = 18, height = 12, units = "cm", res = 200)
varImpPlot(rf.all)
dev.off()


map <- sapply(mks, function(x) {
  mk <- strsplit(x, "\\+")[[1]]
  r1 <- am %in% mk
  names(r1) <- am
  
  rf <- randomForest(fu ~ ., df[names(fu), mk])
  pred <- predict(rf, df.nci[, mk])
  
  identical(names(pred), names(fu2))
  
  r2 <- sapply(unique(dat$nci.xref$too.short), function(too) {
    cor(pred[dat$nci.xref$too.short == too], fu2[dat$nci.xref$too.short == too], 
        use = "pair", method = "spearman")
  })
  
  r3 <- cor(pred, fu2, use= "pair", method = "spearman")
  c(r1, r2, all.too=r3)
})

map <- rbind(map, Roob = oobr)
write.table(round(map, digits = 3), 
            file = "Res/20170317_Fluorouracil/fluTab.absScale.txt", col.names = TRUE, 
            row.names = TRUE, sep = "\t", quote = FALSE)

## Some examples
library(beeswarm)

bdf <- as.data.frame(t(map[16:25, ]))
cc <- unique(dat$nci.xref$color)
names(cc) <- unique(dat$nci.xref$too.short)
cc <- cc[names(bdf)]
cc[length(cc)] <- "pink"

png("Res/20170317_Fluorouracil/absScale/cor.coef.fu.too.tiff", width = 16, height = 12, units = "cm", res = 200)
beeswarm(bdf, col = cc, pch = 19, corral = "gutter", ylab = "correlation coefficent (R)")
bxplot(bdf, add = TRUE)
dev.off()


#########################
sig <- "TYMP+RRM1+UPP1+UCK2+UCK1"
mk <- strsplit(sig, "\\+")[[1]]
rf <- randomForest(fu ~ ., df[names(fu), mk])
pred <- predict(rf, df.nci[, mk])
identical(names(pred), names(fu2))

png("Res/20170317_Fluorouracil/absScale/fu.co.predvsmeasured.3.tiff", 
    width = 20, height = 10, units = "cm", res = 200)
layout(matrix(1:2, 1, 2))
plot(pred[dat$nci.xref$too.short == "CO"], fu2[dat$nci.xref$too.short == "CO"], 
     col = "orange3", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured GI50 (uM)", main = sig)
text(pred[dat$nci.xref$too.short == "CO"], fu2[dat$nci.xref$too.short == "CO"], 
     gsub("_NCI60", "", names(pred[dat$nci.xref$too.short == "CO"])), cex = 0.7)
abline(a = 0, b = 1)

plot(pred[dat$nci.xref$too.short == "LE"], fu2[dat$nci.xref$too.short == "LE"], 
     col = "cyan", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured GI50 (uM)", main = sig)
abline(a = 0, b = 1)
text(pred[dat$nci.xref$too.short == "LE"], fu2[dat$nci.xref$too.short == "LE"], 
     gsub("_NCI60", "", names(pred[dat$nci.xref$too.short == "LE"])), cex = 0.7)
dev.off()

#########################
sig <- "TYMP+RRM1+PPAT+UPP1+UCK1"
mk <- strsplit(sig, "\\+")[[1]]
rf <- randomForest(fu ~ ., df[names(fu), mk])
pred <- predict(rf, df.nci[, mk])
identical(names(pred), names(fu2))

png("Res/20170317_Fluorouracil/absScale/fu.co.predvsmeasured.3.2.tiff", 
    width = 20, height = 10, units = "cm", res = 200)
layout(matrix(1:2, 1, 2))
plot(pred[dat$nci.xref$too.short == "CO"], fu2[dat$nci.xref$too.short == "CO"], 
     col = "orange3", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured GI50 (uM)", main = sig)
text(pred[dat$nci.xref$too.short == "CO"], fu2[dat$nci.xref$too.short == "CO"], 
     gsub("_NCI60", "", names(pred[dat$nci.xref$too.short == "CO"])), cex = 0.7)
abline(a = 0, b = 1)

plot(pred[dat$nci.xref$too.short == "LE"], fu2[dat$nci.xref$too.short == "LE"], 
     col = "cyan", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured GI50 (uM)", main = sig)
text(pred[dat$nci.xref$too.short == "LE"], fu2[dat$nci.xref$too.short == "LE"], 
     gsub("_NCI60", "", names(pred[dat$nci.xref$too.short == "LE"])), cex = 0.7)
abline(a = 0, b = 1)
dev.off()

#########################
sig <- "TYMS+UMPS.1+TYMP"
mk <- strsplit(sig, "\\+")[[1]]
rf <- randomForest(fu ~ ., df[names(fu), mk])
pred <- predict(rf, df.nci[, mk])
identical(names(pred), names(fu2))

png("Res/20170317_Fluorouracil/absScale/fu.co.predvsmeasured.5.tiff", 
    width = 30, height = 10, units = "cm", res = 200)
layout(matrix(1:3, 1, 3))
plot(pred[dat$nci.xref$too.short == "CO"], fu2[dat$nci.xref$too.short == "CO"], 
     col = "orange3", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured GI50", main = sig)
abline(a = 0, b = 1)
plot(pred[dat$nci.xref$too.short == "LE"], fu2[dat$nci.xref$too.short == "LE"], 
     col = "cyan", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured GI50 ((uM)", main = sig)
abline(a = 0, b = 1)
plot(pred[dat$nci.xref$too.short == "BR"], fu2[dat$nci.xref$too.short == "BR"], 
     col = "blue", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured GI50 (uM)", main = sig)
abline(a = 0, b = 1)

dev.off()

#######################
sig <- "TYMS+UMPS.1+UPP1"
mk <- strsplit(sig, "\\+")[[1]]
rf <- randomForest(fu ~ ., df[names(fu), mk])
pred <- predict(rf, df.nci[, mk])
identical(names(pred), names(fu2))

png("Res/20170317_Fluorouracil/absScale/fu.co.predvsmeasured.5.2.tiff", 
    width = 30, height = 10, units = "cm", res = 200)

layout(matrix(1:3, 1, 3))
plot(pred[dat$nci.xref$too.short == "CO"], fu2[dat$nci.xref$too.short == "CO"], 
     col = "orange3", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured GI50 (uM)", main = sig)
abline(a = 0, b = 1)

plot(pred[dat$nci.xref$too.short == "LE"], fu2[dat$nci.xref$too.short == "LE"], 
     col = "cyan", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "mmeasured GI50 (uM)", main = sig)
abline(a = 0, b = 1)

plot(pred[dat$nci.xref$too.short == "BR"], fu2[dat$nci.xref$too.short == "BR"], 
     col = "blue", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured GI50 (uM)", main = sig)
abline(a = 0, b = 1)
dev.off()

####
png("Res/20170317_Fluorouracil/absScale/nci60.gi50.png", width = 28, height = 9, 
    res = 200, units = "cm")
ord <- order(fu2)
layout(matrix(c(1, 1, 2), 1, 3))
barplot(fu2[ord], col = dat$nci.xref$color[ord], las = 2, ylab = "GI50 (uM)", 
          names.arg = gsub("_NCI60", "", names(fu2[ord])))
beeswarm(fu2 ~ dat$nci.xref$too.short, col = cc[levels(factor(dat$nci.xref$too.short))], pch = 19, 
         corral = "gutter", ylab = "GI50 (uM)")
dev.off()


png("Res/20170317_Fluorouracil/absScale/crc65.gi50.png", width = 18, height = 9, 
    res = 200, units = "cm")
ord <- order(fu)
barplot(fu[ord], col = "orange3", las = 2, ylab = "GI50 (uM)", 
        names.arg = gsub("_NCI60", "", names(fu2[ord])))
dev.off()



map2 <- sapply(mks, function(x) {
  mk <- strsplit(x, "\\+")[[1]]
  r1 <- am %in% mk
  names(r1) <- am
  
  rf <- randomForest(fu ~ ., df[names(fu), mk])
  pred <- predict(rf, df.nci[, mk])
  
  identical(names(pred), names(fu2))
  
  r2 <- sapply(unique(dat$nci.xref$too.short), function(too) {
    cor(pred[dat$nci.xref$too.short == too], fu2[dat$nci.xref$too.short == too], 
        use = "pair", method = "pearson")
  })
  
  r3 <- cor(pred, fu2, use= "pair", method = "pearson")
  c(r1, r2, all.too=r3)
})