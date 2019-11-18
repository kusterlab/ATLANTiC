library(TCGAbiolinks)
library(matrixStats)


clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
clinical2 <- GDCquery_clinic(project = "TCGA-READ", type = "clinical")


rf <- function(x) {
  a <- apply(cbind(x$days_to_last_follow_up, x$days_to_death), 2, as.numeric)
  sd <- rowMaxs(a, na.rm = TRUE)
  sd[is.infinite(sd)] <- NA
  data.frame(status = x$vital_status, 
             time = sd, 
             stage = x$tumor_stage,
             disease = x$disease,
             row.names = as.character(x$bcr_patient_barcode))
}
clidat <- rbind(rf(clinical), rf(clinical2))
saveRDS(clidat, file = "Res/20170317_Fluorouracil/clidat.cptac.colon.RDS")

### drug info

query1 <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical")
GDCdownload(query1)
clinical <- GDCprepare_clinic(query, clinical.info = "")