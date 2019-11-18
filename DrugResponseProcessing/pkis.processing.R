library(reshape2)

pkiDrugResponse <- list()
pkiDrugResponse <- within(pkiDrugResponse, {
  pkis <- read.delim("Dat/ElkinsPKIs/nci60_pkis.txt")
  array <- acast(pkis, X ~ CELL ~ ENDPOINT, value.var = "NLOGVALUE", fun.aggregate = mean)
  gi50 <- array[, , "GI50"]
  ic50 <- array[, , "IC50"]
  lc50 <- array[, , "LC50"]
  tgi <- array[, , "TGI"]
  rm(list=c("pkis", "array"))
})



# summary(pkiDrugResponse)
# 
# rownames(pkiDrugResponse$lc50)
# colnames(pkiDrugResponse$lc50)
# head(pkiDrugResponse$annot)


