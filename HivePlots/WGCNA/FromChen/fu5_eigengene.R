library(ROCR)

tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")

k <- readRDS('fromMartin/MElightgoldenrodyellow.rds')
k <- structure(k$MElightgoldenrodyellow, names = rownames(k))
outpath <- "/media/kusterlab/users_files/Martin Frejno/postdoc/projects/phosphoproject/res/20180302_wgcna_allsites_allproteins_incnas_crc65_nonas_crc65/"
ME_crc65_noimputation <- readRDS(file.path(outpath, 'RDS', 'mergedMEs.rds'))$NoImputation
rownames(ME_crc65_noimputation) <- names(k)

msi <- structure(tryp$crc65_xref$MSI.stable, names = tryp$crc65_xref$name_MQ)


####
crc65.resp <- readRDS("fromMartin/5FU.rda")
nn <- paste(rownames(crc65.resp), "_CRC65", sep = "")
nn[nn == "HCT 116_CRC65"] <- "HCT116_CRC65"
nn[nn == "HT-29_CRC65"] <- "HT29_CRC65"
crc65.resp <- structure(crc65.resp[, 1], names = nn)

sens1 <- crc65.resp <= crc65.resp["HRA-19_CRC65"]
sens2 <- crc65.resp <= 5
sens3 <- crc65.resp <= crc65.resp["SK-CO-1_CRC65"]

# par(mar = c(8, 4, 1, 1))
# barplot(sort(crc65.resp), las = 2)

14500917/10000^2


ii <- intersect(names(k), names(msi))

df <- data.frame(sens1 = sens1[ii], 
                 sens2 = sens2[ii],
                 sens3 = sens3[ii],
                 k = k[ii], 
                 msi = msi[ii])
df <- df[!is.na(df$sens1), ]
var <- ME_crc65_noimputation[rownames(df), ]

# ii2 <- crc65.resp < crc65.resp["HRA-19_CRC65"] | crc65.resp > crc65.resp["SK-CO-1_CRC65"]
# df <- df[ii2, ]
# var <- var[ii2, ]

###
callauc <- function(v1, plot = FALSE) {
  
  pred <- prediction(cbind(v1, v1, v1), df[, c("sens1", "sens2", "sens3")])
  perf <- performance(pred, measure="tpr", x.measure="fpr")
  AUC <- as.numeric(performance(pred, measure="auc")@y.values) 
  
  if (plot) {
    
    pred1<-prediction( v1,  df$msi)
    perf1 <- performance(pred1, measure="tpr", x.measure="fpr")
    plot(perf1, xlim = c(0, 1), ylim = c(0, 1), xlab = "False positive rate", ylab = "True positive rate")
    
    num <- as.numeric(df$msi)
    pred0<-prediction(cbind(num, num, num), df[, c("sens1", "sens2", "sens3")])
    perf0 <- performance(pred0, measure="tpr", x.measure="fpr")
    
    for (i in 1:3) {
      lines(perf@x.values[[i]], perf@y.values[[i]], col = i+1, lty = 1)
      points(perf0@x.values[[i]][2], perf0@y.values[[i]][2], pch = 19, col = i+1)
    }
    legend("bottomright", col = 1:4, lty = 1, bty = "n",
           legend = c("var vs MSI", "var vs Bodmer sens", "var vs <5", "var vs Bodmer sens+inter"))  
  }
  AUC
}


callauc(df$k, plot = TRUE)

aucs <- sapply(var, callauc)

barplot(aucs[1, ])
max(aucs[1, ])
callauc(var[, which.max(aucs[1, ])], plot = TRUE)

barplot(aucs[2, ])
max(aucs[2, ])
callauc(var[, which.max(aucs[2, ])], plot = TRUE)
which.max(aucs[2, ])

barplot(aucs[3, ])
max(aucs[3, ])
callauc(var[, which.max(aucs[3, ])], plot = TRUE)


barplot(aucs[, "violetred3"])
max(aucs[1, ])
callauc(var$MEvioletred4, plot = TRUE)
