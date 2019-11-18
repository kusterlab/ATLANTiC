library(openxlsx)
library(TCGAbiolinks)
library(survHD)


###
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")

ii <- grep("fluo|Flour|FU", levels(clinical.drug$drug_name), ignore.case = TRUE, value = TRUE)[-c(4, 5, 7)]

ftreat <- clinical.drug$drug_name %in% ii
drug.fu <- clinical.drug[ftreat, ]


##
query2 <- GDCquery(project = "TCGA-STAD", 
                   data.category = "Clinical")
GDCdownload(query2)
clinical2 <- GDCprepare_clinic(query2, clinical.info = "patient")
clinical2.drug <- GDCprepare_clinic(query2, clinical.info = "drug")

grep("Flour|FU", levels(clinical2.drug$drug_name), ignore.case = TRUE, value = TRUE)

ftreat <- clinical2.drug$drug_name %in% c("5-Flourouracil", "5-FU", "5FU","5-FU+ etoposidium", "FLOUROURACIL")
drug.fu <- clinical2.drug[ftreat, ]

###
##
query3 <- GDCquery(project = "TCGA-ESCA", 
                   data.category = "Clinical")
GDCdownload(query3)
clinical3 <- GDCprepare_clinic(query3, clinical.info = "patient")
clinical3.drug <- GDCprepare_clinic(query3, clinical.info = "drug")

grep("Flour|FU", levels(clinical3.drug$drug_name), ignore.case = TRUE, value = TRUE)

ftreat <- clinical3.drug$drug_name %in% c("5-Flourouracil", "5 FU", "5-FU", "5FU","5-FU+ etoposidium", "FLOUROURACIL")
drug.fu <- clinical3.drug[ftreat, ]


##
query4 <- GDCquery(project = "TCGA-HNSC", 
                   data.category = "Clinical")
GDCdownload(query4)
clinical4 <- GDCprepare_clinic(query4, clinical.info = "patient")
clinical4.drug <- GDCprepare_clinic(query4, clinical.info = "drug")

grep("Flour|FU", levels(clinical4.drug$drug_name), ignore.case = TRUE, value = TRUE)

ftreat <- clinical4.drug$drug_name %in% c("5-Flourouracil", "5 FU", "5-FU", "5FU","5-FU+ etoposidium", "FLOUROURACIL")
drug.fu <- clinical4.drug[ftreat, ]


###
query5 <- GDCquery(project = "TCGA-COAD", 
                   data.category = "Clinical")
GDCdownload(query5)
clinical5 <- GDCprepare_clinic(query5, clinical.info = "patient")
clinical5.drug <- GDCprepare_clinic(query5, clinical.info = "drug")

grep("Flour|FU", levels(clinical5.drug$drug_name), ignore.case = TRUE, value = TRUE)

ftreat <- clinical5.drug$drug_name %in% c("5-Flourouracil", "5- FU", "oxaliplatinum+ 5-FU", "5 FU", "5-FU", "5FU","5-FU+ etoposidium", "FLOUROURACIL")
drug.fu <- clinical5.drug[ftreat, ]


###
query6 <- GDCquery(project = "TCGA-READ", 
                   data.category = "Clinical")
GDCdownload(query6)
clinical6 <- GDCprepare_clinic(query6, clinical.info = "patient")
clinical6.drug <- GDCprepare_clinic(query6, clinical.info = "drug")



ftreat <- clinical6.drug$drug_name %in% grep("Fluor|FU", levels(clinical6.drug$drug_name), ignore.case = TRUE, value = TRUE)
drug.fu <- clinical6.drug[ftreat, ]

clinical.drug$CancerType <- "BRCA"
clinical2.drug$CancerType <- "STAD"
clinical3.drug$CancerType <- "ESCA"
clinical4.drug$CancerType <- "HNSC"
clinical5.drug$CancerType <- "COAD"
clinical6.drug$CancerType <- "READ"
cli.drug <- rbind(clinical.drug, clinical2.drug, clinical3.drug, 
                  clinical4.drug, clinical5.drug, clinical6.drug)

ii <- grep("fluo|flour|fu", levels(cli.drug$drug_name), ignore.case = TRUE, value = TRUE)
ii <- ii[!grepl("Ful", ii, ignore.case = TRUE)]
fu.drug <- cli.drug[cli.drug$drug_name %in% ii, ]

fu5 <- list(all6cancer = cli.drug, 
            fu5 = fu.drug)

saveRDS(fu5, file = "Res/20170324_Fluororacil_TCGA/cli.drug.RDS")
