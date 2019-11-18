library(stringr)
library(stringi)

# fits <- c("crc65.ccle.allsites",
#           "nci60.ccle.regsites.gluc",
#           "nci60.dtp.regsites.gluc",
#           "crc65.ccle.proteinGroups",
#           "nci60.ccle.regsites",
#           "nci60.dtp.regsites",
#           "crc65.ctrp.proteinGroups",
#           "nci60.ctrp.allsites.gluc",
#           "nci60.gdsc.allsites.gluc",
#           "crc65.gdsc.allsites",
#           "nci60.ctrp.allsites",
#           "nci60.gdsc.allsites",
#           "crc65.gdsc.proteinGroups",
#           "nci60.ctrp.proteinGroups",
#           "nci60.gdsc.proteinGroups",
#           "crc65.medico.proteinGroups",
#           "nci60.ctrp.regsites.gluc",
#           "nci60.gdsc.regsites.gluc",
#           "nci60.ctrp.regsites",
#           "nci60.gdsc.regsites",
#           "nci60.ccle.allsites.gluc",
#           "nci60.dtp.allsites.gluc",
#           "nci60.ccle.allsites",
#           "nci60.dtp.allsites",
#           "nci60.ccle.proteinGroups",
#           "nci60.dtp.proteinGroups")

fits <- c("crc65.ccle.allsites",
          "crc65.ccle.proteinGroups",
          "crc65.ccle.regsites",
          "crc65.ctrp.allsites",
          "crc65.ctrp.proteinGroups",
          "crc65.ctrp.regsites",
          "crc65.gdsc.allsites",
          "crc65.gdsc.proteinGroups",
          "crc65.gdsc.regsites",
          "crc65.medico.allsites",
          "crc65.medico.proteinGroups",
          "crc65.medico.regsites",
          "nci60.ccle.allsites.gluc",
          "nci60.ccle.allsites",
          "nci60.ccle.proteinGroups",
          "nci60.ccle.regsites.gluc",
          "nci60.ccle.regsites",
          "nci60.ctrp.allsites.gluc",
          "nci60.ctrp.allsites",
          "nci60.ctrp.proteinGroups",
          "nci60.ctrp.regsites.gluc",
          "nci60.ctrp.regsites",
          "nci60.dtp.allsites.gluc",
          "nci60.dtp.allsites",
          "nci60.dtp.proteinGroups",
          "nci60.dtp.regsites.gluc",
          "nci60.dtp.regsites",
          "nci60.gdsc.allsites.gluc",
          "nci60.gdsc.allsites",
          "nci60.gdsc.proteinGroups",
          "nci60.gdsc.regsites.gluc",
          "nci60.gdsc.regsites")

elnetlist <- list()
for (i in fits) {
  print(i)
  fit <- readRDS(paste("Res/20170110_mq15processLOD2/", i, ".RDS", sep = ""))
  a <- lapply(names(fit$fit), function(x) {
    xx <- fit$fit[[x]]
    if (class(xx) == "try-error")
      return()
    if (is.null(xx[1]))
      return()
    if (is.na(xx[1]))
      return()
    sum <- as.data.frame(xx$summary)
    sum$drug <- paste(x, i)
    selectCol <- c("drug", "gene", "correlation", "meanefsize", "sdefsize", "nonzerocounts")
    sum[, selectCol]
  })
  names(a) <- names(fit$fit)
  elnetlist[[i]] <- a[sapply(a, length) != 0]
}
sapply(elnetlist, length)

##
allres <- unlist(elnetlist, recursive = FALSE)
allres <- lapply(allres, function(x) {
  x <- x[x$nonzerocounts > 6, c("drug", "gene", "meanefsize", "nonzerocounts")]
  x$meanefsize <- scale(x$meanefsize, center = FALSE, scale = TRUE)
  x
})
allres <- do.call("rbind", allres)

head(allres$drug)

##
rss <- str_split(stri_reverse(allres$drug), " ", n = 2, simplify = TRUE)
d <- rss <- apply(rss, 2, stri_reverse)
allres$drug <- d[, 2]
allres$inputData <- d[, 1]

igluc <- grep(".gluc", allres$inputData)
allres$inputData[-igluc] <- paste(allres$inputData[-igluc], "trypsin", sep=".")

allres$inputData <- gsub(".ccle.|.ctrp.|.gdsc.|.medico.|.dtp.", ".", allres$inputData)
rownames(allres) <- NULL

object.size(allres)
# table(allres$inputData)


head(allres)

r2 <- stri_reverse(allres$drug)
r2 <- str_split(r2, pattern = "_", n=2, simplify = TRUE)
r2 <- stri_reverse(r2[, 1])
table(r2)

unique(allres$inputData)

allres$inputData2 <- NULL

##
ff <- factor(allres$inputData2, 
             labels = c("CRC65.PSITE.TRYP.ALL", 
                        "CRC65.PROT.TRYP", 
                        "CRC65.PSITE.TRYP.REGSITE",
                        "NCI60.PSITE.GLUC.ALL", 
                        "NCI60.PSITE.TRYP.ALL",
                        "NCI60.PROT.TRYP",
                        "NCI60.PSITE.GLUC.REGSITE", 
                        "NCI60.PSITE.TRYP.REGSITE"))




rr <- data.frame(
  EFFECTSIZE = allres$meanefsize, 
  FREQ = allres$nonzerocounts,
  ID = allres$gene, 
  DRUG = allres$drug, 
  DRUGDATA = r2, 
  INPUTDATA = ff)
rr$ID <- as.character(rr$ID)
rr$DRUG <- as.character(rr$DRUG)


levels(rr$DRUGDATA)

# rr <- readRDS(file = "Res/20170112_elnetResSummary/elnetsSummary.RDS")

tryp <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170110_mq15processLOD2/sites.trypsin.RDS")
gluc <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Res/20170110_mq15processLOD2/sites.gluC.RDS")

rr$AA <- NA
rr$POSITION <- NA
for (i in c(".PSITE.TRYP", ".PSITE.GLUC")) {
  idx <- grepl(i, rr$INPUTDATA)
  annot <- switch (i,
                   ".PSITE.TRYP" = tryp$annot,
                   ".PSITE.GLUC" = gluc$annot)
  aa <-substr(annot$`Sequence window`, 16, 16)
  names(aa) <- rownames(annot)
  rr$AA[idx] <- aa[as.character(rr$ID)[idx]]
  
  bb <- annot$Position
  names(bb) <- rownames(annot)
  rr$POSITION[idx] <- bb[as.character(rr$ID)[idx]]
}

head(rr)

saveRDS(rr, file = "Res/20170112_elnetResSummary/elnetsSummary.RDS")



##

aa <- allres[allres$drug == "743414_Imatinib_DTP" & allres$inputData == "nci60.regsites.trypsin", ]
plot(aa$meanefsize, aa$nonzerocounts)
aa <- allres[allres$drug == "imatinib" & allres$inputData == "nci60.regsites", ]
plot(aa$meanefsize, aa$nonzerocounts)




