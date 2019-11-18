library(stringr)
library(fastmatch)

kin <- read.delim("Dat/Uniprot/Annotation_Kinase_list.txt", stringsAsFactors = FALSE)

tryp <- readRDS('/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/dat/preproc/crc65_psites_annotation.rds')
fp <- readRDS('/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/dat/preproc/crc65_fp_annotation.rds')

mod_crc <- read.csv("/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/res/20180302_wgcna_allsites_allproteins_incnas_crc65_nonas_crc65/CSVs/wgcna_NoImputation.csv", stringsAsFactors = FALSE)
mod_nci <- read.csv("/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/res/20180308_wgcna_allsites_allproteins_incnas_nci60_nonas_nci60/CSVs/wgcna_NoImputation.csv", stringsAsFactors = FALSE)


# col - "Gene", "id", "sequence window", "position", "isKinase", "panel"
mod_crc$gene <- str_split_fixed(mod_crc$id, "_", 2)[, 1]
mod_crc$isKinase <- ""
mod_crc$isKinase[mod_crc$gene %in% kin$Name] <- "TRUE"
mod_crc$sequenceWindow <- tryp$`Sequence window`[fmatch(mod_crc$id, tryp$label)]
mod_crc$position <- tryp$Positions[fmatch(mod_crc$id, tryp$label)]
mod_crc$panel <- "CRC65"


mod_nci$gene <- str_split_fixed(mod_nci$id, "_", 2)[, 1]
mod_nci$isKinase <- ""
mod_nci$isKinase[mod_nci$gene %in% kin$Name] <- "TRUE"
mod_nci$sequenceWindow <- tryp$`Sequence window`[fmatch(mod_nci$id, tryp$label)]
mod_nci$position <- tryp$Positions[fmatch(mod_nci$id, tryp$label)]
mod_nci$panel <- "NCI60"

mod <- rbind(mod_nci, mod_crc)
mod <- mod[, c("gene", "isKinase", "id", "moduleColor", "sequenceWindow", "position", "panel")]

head(mod)
write.table(mod, file = "Res/20180314_wgcnaSupTables/moduleCandidates.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

