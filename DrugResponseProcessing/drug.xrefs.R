## Total number study
# both NCI60 and CRC65 - CTRP, GDSC, CCLE
# NCI60 only - DTP, PKIS
# CRC only - cetuximab

###############################################################
nci60 <- readRDS("../martin/Drug_sensitivity/drugs.rda")
head(nci60)
crc <-readRDS("../martin/Drug_sensitivity/drugs_crc.rda")
head(crc)
spl <- strsplit(as.character(crc$drug), split = "_")
dd <- do.call(rbind, spl)
head(dd)
crc$dataset <- dd[, ncol(dd)]
crc$drug <- gsub("_CCLE|_CTRP|_GDSC|_MEDICO", "", crc$drug)
head(crc)
dtp <- read.delim("~/CloudChen/Projects/PhosphoNCI60/Data/DTP_DrugSensitivity/DTP_NCI60.xls/DTP_NCI60_1_5_1.txt", 
                  comment.char = "#", check.names = FALSE, as.is = TRUE)
kbs <- read.delim("Dat/kbscreen&review/kinobeadScreen.txt")
santosRev <- read.delim("Dat/kbscreen&review/santos_nrd.2016.230-s2.txt")

############################## combine list #################################
drug <- list(santosRev = santosRev$PARENT_PREF_NAME, 
             kinoScreen = kbs$Drug,
             dtp = dtp$Name, 
             dtp.auc = as.character(nci60$drug[nci60$dataset=="DTP"]),
             CCLE = as.character(nci60$drug[nci60$dataset=="CCLE"]),
             CTRP = as.character(nci60$drug[nci60$dataset=="CTRP"]),
             GDSC = as.character(nci60$drug[nci60$dataset=="GDSC"]))
drug <- lapply(drug, as.character)
drug <- lapply(drug, unique)
drug.l <- lapply(drug, tolower)
drug.nop <- lapply(drug.l, function(x) gsub("[[:punct:]]| ", "", x))

############################## for santos review #################################
srmat <- sapply(drug.nop[c(1, 3:7)], function(x) drug.nop$santosRev %in% x)
rownames(srmat) <- drug.nop$santosRev
srmat <- srmat[rowSums(srmat) > 1, ]
colSums(srmat)

res <- matrix(NA, nrow(srmat), ncol(srmat), 
              dimnames = list(rownames(srmat), colnames(srmat)))
cn <- colnames(srmat)
for (i in cn) {
  si <- drug.nop[[i]] %in% rownames(srmat)[srmat[, i]]
  res[srmat[, i], i] <- sort(drug[[i]][si])
}

apply(res, 2, function(x) {
  ir <- !is.na(x)
  all(rownames(res)[ir] == gsub("[[:punct:]| ]", "", tolower(x[ir])))
})

write.table(res, file = "Dat/xref/xref.rlq.santosRev.txt", 
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)



############################## for kinobead drug screen #################################
srmat <- sapply(drug.nop[c(2, 3:7)], function(x) drug.nop$kinoScreen %in% x)
rownames(srmat) <- drug.nop$kinoScreen
srmat <- srmat[rowSums(srmat) > 1, ]
colSums(srmat)

res <- matrix(NA, nrow(srmat), ncol(srmat), 
              dimnames = list(rownames(srmat), colnames(srmat)))
cn <- colnames(srmat)
for (i in cn) {
  si <- drug.nop[[i]] %in% rownames(srmat)[srmat[, i]]
  res[srmat[, i], i] <- sort(drug[[i]][si])
}

apply(res, 2, function(x) {
  ir <- !is.na(x)
  all(rownames(res)[ir] == gsub("[[:punct:]]| ", "", tolower(x[ir])))
})

write.table(res, file = "Dat/xref/xref.rlq.kinoScreen.txt", 
            col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)












