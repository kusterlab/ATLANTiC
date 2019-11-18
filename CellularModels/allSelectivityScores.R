library(data.table)
library(parallel)
source("functionsForShiny.R")
dd <- readRDS("../atlantic_old/celllinatorshiny/dat_directBinders.RDS")


at <- lapply(names(dd)[1:4], function(nam) {
  
  curveinformation <- NULL
  value <- list()
  value$drgTotal <- dd$cpd
  desired_drug_total <- c(colnames(value$drgTotal[-c(1:3)]))
  
  value$exprTotal <- list(inputData1 = dd[[nam]])
  rs <- mclapply(value$drgTotal$Gene.names, mc.cores = 20, function(desired_protein) {
    all_protein_drug_ss_list  <- lapply(desired_drug_total, function(desired_drug , expression_list, drug, desired_protein) {
      drug_protein_ss_calc(desired_drug =  desired_drug,
                           expression_list = expression_list , 
                           drug =  drug , 
                           desired_protein =  desired_protein, 
                           curveinformation = curveinformation)
    } , expression_list  = value$exprTotal, drug = value$drgTotal, desired_protein =  desired_protein)
    as.data.frame(rbindlist(all_protein_drug_ss_list), stringsAsFactors = F)
  })
  rs <- do.call(rbind, rs)
  rs$Dataset <- nam
  rs
}) 

at1 <- do.call(rbind, at)
at1 <- at1[order(at1$Selectivity_Score, decreasing = TRUE), ]
for (i in colnames(at1)) {
  if (is.character(at1[[i]]))
    at1[[i]] <- factor(at1[[i]])
}
at1$Selectivity_Score <- signif(at1$Selectivity_Score, digits = 3)
saveRDS(at1, file = "../../cellinator2/allSelectivityScores.RDS")

###
at1 <- readRDS("../../cellinator2/allSelectivityScores.RDS")

dd$selectivityScore <- at1
###
value <- list()
value$drgTotal <- dd$cpd
desired_drug_total <- c(colnames(value$drgTotal[-c(1:3)]))
value$drgTotal$Gene.name


dres <- dd[1:4]
for (i in 1:length(dres)) {
  dres[[i]] <- z_scored(dres[[i]])
  dres[[i]] <- dres[[i]][dres[[i]]$Gene.names %in% value$drgTotal$Gene.names, ]
  dres[[i]] <- setDT(melting(dres[[i]]))
  dres[[i]]$Dataset <- names(dres)[i]
}
dd$zscore <- dres

saveRDS(dd, file = "../../cellinator2/dat.RDS")

####
library(data.table)
dd <- readRDS(file = "../../cellinator2/dat.RDS")
for (i in 1:length(dd$zscore)) {
  dd$zscore[[i]] <- setDT(dd$zscore[[i]])
}
saveRDS(dd, file = "../../cellinator2/dat.RDS")

