library(reshape2)

drug <- readRDS("../martin/Drug_sensitivity/dugs_crc_log.rda")
drug$dataset <- sapply(strsplit(as.character(drug$drug), "_"), function(x) x[length(x)])

extractData <- function(dataset, drug) {
  i <- drug$parameter == "AUC" & drug$dataset == dataset
  auc <- drug[ i, ]
  assay <- acast(auc, cellline ~ drug, value.var = "value")
  t(assay)
}

ds <- c("CTRP", "GDSC", "CCLE", "MEDICO")
aucs <- lapply(ds, extractData, drug)
names(aucs) <- ds


slist <- readRDS("Res/20161230_mq15finalprocess/sites.RDS")

#
names <- gsub("_CRC65", "", colnames(slist$crc65_expr))
fitCols <- function(drugm, names) {
  mat <- matrix(NA, nrow = nrow(drugm), ncol = length(names), 
                dimnames = list(rownames(drugm), names))
  ii <- intersect(names, colnames(drugm))
  mat[, ii] <- drugm[, ii]
  1-mat
}

aucsfit <- lapply(aucs, fitCols, names = names)
drug.crc65 <- aucsfit

rm(list = c("ds", "drug", "extractData", "slist", 'names', "fitCols", "aucsfit", "aucs"))
# sapply(aucs, dim)


