library(reshape2)
drug <- read.csv("../martin/Drug_sensitivity/drugs.csv")
head(drug)

auc <- drug[ drug$parameter == "AUC", ]

assay <- acast(auc, cellline ~ drug ~ dataset, value.var = "value")
dim(assay)
levels(auc$dataset)

ccle <- assay[, , "CCLE"]
ctrp <- assay[, , "CTRP"]
gdsc <- assay[, , "GDSC"]

ccle <- t(ccle[, colSums(!is.na(ccle)) >= 24]) # lapatinib
ctrp <- t(ctrp[, colSums(!is.na(ctrp)) >= 24]) # has imatinib, lapatinib
gdsc <- t(gdsc[, colSums(!is.na(gdsc)) >= 24]) # 

rm(list=c("auc", "assay", "drug"))

# grep("lapatinib", colnames(ccle), ignore.case = TRUE, value = TRUE)
# grep("lapatinib", colnames(ctrp), ignore.case = TRUE, value = TRUE)
# grep("lapa", colnames(gdsc), ignore.case = TRUE, value = TRUE)
# 
# 
# write.table(colnames(ccle), file = "temp1.txt", col.names = TRUE, row.names = TRUE, quote=FALSE, sep="\t")
# write.table(drug$annot.dtp, file = "temp2.txt", col.names = TRUE, row.names = TRUE, quote=FALSE, sep="\t")



