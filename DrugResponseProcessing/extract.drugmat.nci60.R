xref <- read.delim("Dat/xref/xref.cellline.nci60.txt", as.is = TRUE)
source("R/drugResponseProcess/aucs.selection&process.R")
source("R/drugResponseProcess/dtp.selection&process.R")

mapNCI60 <- function(drug.mat, map.name, res.name) {
  nr <- length(map.name)
  
  m <- matrix(NA, nrow(drug.mat), nr, dimnames = list(rownames(drug.mat), map.name))
  ina <- intersect(colnames(drug.mat), map.name)
  m[, ina] <- drug.mat[, ina]
  colnames(m) <- res.name
  m
}

drug.nci60 <- list()
drug.nci60$dtp.annot <- dtp$dtp.annot
drug.nci60$dtp.gi50 <- mapNCI60(dtp$dtp.gi50, xref$name_DTP, xref$name_MQ)
drug.nci60$ccle.auc <- 1 - mapNCI60(aucs$CCLE, xref$name_AUC, xref$name_MQ)
drug.nci60$gdsc.auc <- 1 - mapNCI60(aucs$GDSC, xref$name_AUC, xref$name_MQ)
drug.nci60$ctrp.auc <- 1 - mapNCI60(aucs$CTRP, xref$name_AUC, xref$name_MQ)

rownames(drug.nci60$ccle.auc) <- paste(rownames(drug.nci60$ccle.auc), "CCLE", sep = "_")
rownames(drug.nci60$gdsc.auc) <- paste(rownames(drug.nci60$gdsc.auc), "GDSC", sep = "_")
rownames(drug.nci60$ctrp.auc) <- paste(rownames(drug.nci60$ctrp.auc), "CTRP", sep = "_")
rownames(drug.nci60$dtp.gi50) <- paste(rownames(drug.nci60$dtp.gi50), "DTP", sep = "_")
rownames(drug.nci60$dtp.annot) <- paste(rownames(drug.nci60$dtp.annot), "DTP", sep = "_")

rm(list  = c("xref", "aucs", "dtp", "mapNCI60"))

# 
# rownames(gdsc.auc)
# barplot(sort(gdsc.auc["imatinib", ]), las = 2)
# barplot(sort(ctrp.auc["imatinib", ]), las = 2)
# 
# ii <- intersect(rownames(gdsc.auc), rownames(ctrp.auc))
# 
# for (iii in ii) {
#   plot(gdsc.auc[iii, ], ctrp.auc[iii, ], main = iii)
# }

