dtpMOA <- list()
suppressWarnings(dtpMOA <- within(dtpMOA, {
  dtp <- read.delim("/home/chen/CloudChen/Projects/PhosphoNCI60/Data/DTP_DrugSensitivity/DTP_NCI60.xls/DTP_NCI60_1_5_1.txt", 
                    comment.char = "#")
  # idx <- dtp$Total.after.quality.control >= 2 & !dtp$FDA.status %in% c("", "-")
  idx <- dtp$Total.after.quality.control >= 1 & 
    dtp$MOA %in% c("Ang", "Apo", "CC", "GTPase", "HDAC", "Ho", "MET", "P90", "PI3K", "PKC", "Wnt",
                   "PD", "Pr", "RXR", "STK", "TFG", "Tu", "YK", "ER", "HS", "NFkB", "Pase", "PTK") & 
    !dtp$FDA.status %in% c("", "-")
  idx[is.na(idx)] <- FALSE
  annot.dtp <- dtp[idx, c(1:4, 65:66)]
  rownames(annot.dtp) <- annot.dtp$NSC
  drug.sens <- apply(dtp[idx, c(5:64)], 2, function(x) as.numeric(as.character(x)))
  rownames(drug.sens) <- rownames(annot.dtp)
  remove(list = c("dtp", "idx"))
}))
