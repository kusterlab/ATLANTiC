library(stringr)

kin <- read.delim("Dat/Uniprot/Annotation_Kinase_list.txt", stringsAsFactors = FALSE)

## CRC65
module <- read.csv("/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/res/20180302_wgcna_allsites_allproteins_incnas_crc65_nonas_crc65/CSVs/wgcna_NoImputation.csv", stringsAsFactors = FALSE)
module$kinase <- str_split_fixed(module$id, pattern = "_", 2)[, 1]
module <- module[module$kinase %in% kin$Name, ]
module <- module[c(3, 1, 2)]
module$panel <- "CRC65"
module_crc <- module



## NCI60
module <- read.csv("/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/res/20180308_wgcna_allsites_allproteins_incnas_nci60_nonas_nci60/CSVs/wgcna_NoImputation.csv", stringsAsFactors = FALSE)
module$kinase <- str_split_fixed(module$id, pattern = "_", 2)[, 1]
module <- module[module$kinase %in% kin$Name, ]
module <- module[c(3, 1, 2)]
module$panel <- "NCI60"
module_nci <- module

###
mod <- rbind(module_nci, module_crc)
write.table(mod, file = "Res/20180314_wgcnaSupTables/KinaseModulePresence.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
