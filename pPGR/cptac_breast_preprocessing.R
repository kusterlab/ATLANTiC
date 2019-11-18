greb1 <- read.delim("Dat/cptac/BREAST_GREB1.txt")
v0 <- as.matrix(greb1[, -(1:28)])

x <- v0[1, ]

x <- na.omit(x)
x <- x[grep("TCGA", names(na.omit(x)), value = TRUE)]
names(x) <- paste0("TCGA-", gsub("\\.", "-", substr(grep("TCGA", names(na.omit(x)), value = TRUE), 1, 7)))
x <- x[!duplicated(names(x))]


library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical", 
                  barcode = bc)
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")


identical(names(x), clinical$bcr_patient_barcode)
clinical$GREB1S320 <- x
clinical$OS <- pmax(clinical$days_to_death, clinical$days_to_last_followup, na.rm = TRUE)

dir.create("Res/20180411_breastCptacSurv")
saveRDS(clinical, "Res/20180411_breastCptacSurv/clin.RDS")
saveRDS(clinical.drug, "Res/20180411_breastCptacSurv/clin.drug.RDS")
