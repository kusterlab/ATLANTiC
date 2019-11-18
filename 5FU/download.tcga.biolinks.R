library(RTCGAToolbox)
library(genefilter)
library(TCGAbiolinks)

pat <- readRDS(file = "Res/20170324_Fluororacil_TCGA/tcga.patient.drug.RDS")

bc <- split(pat$bcr_patient_barcode, pat$CancerType)

names(bc)

qv <- lapply(names(bc)[-1], function(i) {
  print(i)
  query <- GDCquery(project = paste0("TCGA-", i),
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq", 
                    file.type  = "normalized_results",
                    experimental.strategy = "RNA-Seq",
                    barcode = bc[[i]],
                    legacy = TRUE)
  GDCdownload(query, method = "api", chunks.per.download = 10)
  GDCprepare(query)
})

