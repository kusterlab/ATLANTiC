library(data.table)
library(fastmatch)
library(stringr)
library(reshape2)


str <- fread("Dat/STRING/9606.protein.links.detailed.v10.5.txt", integer64 = "double", data.table = FALSE)
hist(str$combined_score)
str <- str[str$combined_score > 750, ]
str$protein1 <- gsub("9606.", "", str$protein1)
str$protein2 <- gsub("9606.", "", str$protein2)
str <- str[order(str$combined_score, decreasing = TRUE), ]
idmap <- fread("Dat/Uniprot/HUMAN_9606_idmapping.dat", data.table = FALSE, header = FALSE)
str$uniprotA <- idmap[fmatch(str$protein1, idmap$V3), 1]
str$uniprotB <- idmap[fmatch(str$protein2, idmap$V3), 1]
str$uniprotA <- str_split_fixed(string = str$uniprotA, pattern = "-", n = 2)[, 1]
str$uniprotB <- str_split_fixed(string = str$uniprotB, pattern = "-", n = 2)[, 1]
pt <- paste(idmap$V1, idmap$V2, sep = "_")
str$symbolA <- idmap[fmatch(paste(str$uniprotA, "Gene_Name", sep = "_"), pt), 3]
str$symbolB <- idmap[fmatch(paste(str$uniprotB, "Gene_Name", sep = "_"), pt), 3]

##
mod_crc <- read.csv("/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/res/20180302_wgcna_allsites_allproteins_incnas_crc65_nonas_crc65/CSVs/wgcna_NoImputation.csv", stringsAsFactors = FALSE)
mod_crc$gene <- str_split_fixed(mod_crc$id, "_", 2)[, 1]
mod_crc <- split(mod_crc, mod_crc$moduleColor)
names(mod_crc) <- paste(names(mod_crc), "CRC65")

mod_nci <- read.csv("/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/res/20180308_wgcna_allsites_allproteins_incnas_nci60_nonas_nci60/CSVs/wgcna_NoImputation.csv", stringsAsFactors = FALSE)
mod_nci$gene <- str_split_fixed(mod_nci$id, "_", 2)[, 1]
mod_nci <- split(mod_nci, mod_nci$moduleColor)
names(mod_nci) <- paste(names(mod_nci), "NCI60")

allmod <- c(mod_crc, mod_nci)
ll <- lapply(names(allmod), function(x) {
  mod <- allmod[[x]]s
  s <- str[str$symbolA %in% mod$gene & str$symbolB %in% mod$gene, ]
  if (nrow(s) == 0)
    return(s)
  s$module <- x
  s
})
lt <- do.call(rbind, ll)

md <- str_split_fixed(lt$module, " ", 2)

lt$moduleColor <- md[, 1]
lt$panel <- md[, 2]
lt$module <- NULL
lt <- lt[, c("moduleColor", "panel", "uniprotA", "uniprotB", "symbolA", "symbolB", "protein1", "protein2",
             "neighborhood", "fusion", "cooccurence", "coexpression", "experimental", "database", "textmining", "combined_score")]
head(lt)
write.table(lt, file = "Res/20180314_wgcnaSupTables/moduleString.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
