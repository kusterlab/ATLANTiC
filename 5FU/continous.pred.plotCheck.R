res.lm <- readRDS(file = "Res/20170317_Fluorouracil/linearModelsPredCorTK1phos.RDS")
# boxplot(res.lm)
res.lm <- do.call(cbind, res.lm)
cmr <- colMedians(res.lm)
names(cmr) <- colnames(res.lm)

hist(cmr, breaks = 100)
sort(cmr, decreasing = TRUE)[1:10]

###

res.rf <- readRDS("Res/20170317_Fluorouracil/RFModelsPredCorTK1phos.RDS")

boxplot(res.rf)
which.min(res[[3]])

i <- !is.na(fu)
df <- dat$lod2(dat$nci60.prot)

crc65.resp <- readRDS("fromMartin/5FU.rda")
nn <- paste(rownames(crc65.resp), "_CRC65", sep = "")
nn[nn == "HCT 116_CRC65"] <- "HCT116_CRC65"
nn[nn == "HT-29_CRC65"] <- "HT29_CRC65"
crc65.resp <- structure(crc65.resp[, 1], names = nn)

mks <- names(which.min(res[[5]]))
mks <- strsplit(mks, "\\+")[[1]]
df.crc <- dat$lod2(dat$crc65.prot)

mks <- c("PPAT", "UCK2")