library(reshape2)

drug <- readRDS("../martin/Drug_sensitivity/drugs_log.rda")

extractData <- function(dataset, drug) {
  i <- drug$parameter == "AUC" & drug$dataset == dataset
  auc <- drug[ i, ]
  assay <- acast(auc, cellline ~ drug, value.var = "value")
  assay <- t(assay)
}

# 
# dataset <- unique(drug$dataset)[4]
# table(drug$dataset)
# x <- extractData(dataset = unique(drug$dataset)[1], drug = drug)
# dim(x)
# head(x)

unique(drug$dataset)
ds <- c("CCLE", "CTRP", "GDSC")
aucs <- lapply(ds, extractData, drug)
names(aucs) <- ds

rm(list = c("ds", "drug", "extractData"))






