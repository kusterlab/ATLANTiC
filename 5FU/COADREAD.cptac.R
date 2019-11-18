library(survHD)
library(randomForest)
library(matrixStats)
library(openxlsx)

################################################################################
# LOADING AND PROCESSING DATA
dat <- readRDS(file = "Res/20170317_Fluorouracil/prep.dat.RDS")

# CRC 65
df <- dat$lod2(dat$crc65.prot)
df$ABCC4.1 <- NULL
crc65.resp <- readRDS("fromMartin/5FU.rda")
nn <- paste(rownames(crc65.resp), "_CRC65", sep = "")
nn[nn == "HCT 116_CRC65"] <- "HCT116_CRC65"
nn[nn == "HT-29_CRC65"] <- "HT29_CRC65"
crc65.resp <- structure(crc65.resp[, 1], names = nn)
fu <- crc65.resp[!is.na(crc65.resp)]
barplot(sort(fu), las = 2)

# NCI 60
fu.nci <- (dat$fu19893 + dat$fu757036)/2
df.nci <- dat$lod2(dat$nci60.prot)

# clinical data
cli <- readRDS("Res/20170317_Fluorouracil/clidat.cptac.colon.RDS")



############################# processing CPTAC data ############################
markers <- c("ABCC3", "ABCC4", "TK1", "TYMS", "UMPS", "TYMP", "RRM1", "RRM2", "PPAT", 
             "DPYD", "UPP1", "UCK2")#, "UCK1", "SLC29A1")


cptac <- read.delim("fromMartin/protein_cptac_wo_imputation.txt", check.names = FALSE)
rownames(cptac) <- cptac[, 1]
cptac <- cptac[, -1]
cptac <- as.matrix(cptac)
cptac <- log10(2^cptac)
hist(cptac)
boxplot(cptac)
idmap <- read.csv("/media/msdata5/users_files/Chen/Projects/IntegrativeClustering/Dat/Martin/TCGAmapping.csv")
rownames(idmap) <- idmap$shortID
colnames(cptac) <- substr(as.character(idmap[colnames(cptac), 1]), 1, 12)

cptac <- sweep(cptac, 2, colMedians(cptac, na.rm = TRUE), "-") + median(cptac, na.rm = TRUE)
boxplot(cptac)
##
sub <- dat$lod2(as.data.frame(t(cptac[markers, ])))
boxplot(sub, las = 2)

# df.range <- sapply(df, function(x) max(x) - min(x))
# sub.range <- sapply(sub, function(x) max(x) - min(x))
# fac.range <- df.range[names(sub.range)]/sub.range
# 
# sub <- sweep(sub, 2, fac.range, "*")
sub <- sweep(sub, 2, sapply(df, min)[colnames(sub)], "+")
boxplot(sub, las = 2)


################################################################################
# select a signature/ call from here
################################################################################
# sig <- "TYMS+UMPS+PPAT+DPYD"
sig <- "ABCC3+TYMS+TYMP"

# sig <- "TYMS+RRM1+DPYD+UCK2"



############################## predict on NCI60 data ###########################
mk <- strsplit(sig, "\\+")[[1]]
rf <- randomForest(fu ~ ., df[names(fu), mk])
pred <- predict(rf, df.nci[, mk])
identical(names(pred), names(fu.nci))

plot(pred[dat$nci.xref$too.short == "CO"], fu.nci[dat$nci.xref$too.short == "CO"], 
     col = "orange3", pch = 20, cex  = 2, xlab= "predicted GI50 (uM)", 
     ylab = "measured -log10(GI50)", main = sig)


################################################################################


ftab <- read.delim("Res/20170317_Fluorouracil/fluTab.txt")
mm <- apply(ftab[markers, ], 2, function(x) markers[x == 1])

ss <- sapply(mm, function(mk) {
  rf <- randomForest(fu ~ ., df[names(fu), mk])
  
  pred <- predict(rf, sub[, mk])
  barplot(sort(pred), las = 2)
  names(pred) <- substr(names(pred), 1, 12)
  
  
  survdf <- cli[names(pred), ]
  sur <- Surv(survdf$time, as.integer(survdf$status=="dead"))
  r <- try(plotKMStratifyBy(method = "Lau92", sur, linearriskscore = pred, HRpos = "bottomleft",
                        legendpos = "topright", mark.time = TRUE, censor.at = 1000))
  try(r$hr)
})


unname(t(ss))
####
pred <- predict(rf, sub[, mk])
barplot(sort(pred), las = 2)
names(pred) <- substr(names(pred), 1, 12)


survdf <- cli[names(pred), ]
sur <- Surv(survdf$time, as.integer(survdf$status=="dead"))
r <- plotKMStratifyBy(method = "Lau92", sur, linearriskscore = pred, HRpos = "bottomleft",
                      legendpos = "topright", mark.time = TRUE, censor.at = 1000)
r$hr
r$p.value
