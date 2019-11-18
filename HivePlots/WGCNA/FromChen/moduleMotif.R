library(PTMotif)
library(stringr)

l1 <- readRDS(file = "Res/20180308_motif_with_wo_imputation_crc65/noimptuation_motif_crc65.RDS")
l1 <- do.call(rbind, l1)
l1 <- l1[which(l1$gw.fg.count > 2), ]
l1 <- annotKnownMotifs(l1)
l1$panel <- "CRC65"

l2 <- readRDS(file = "Res/20180312_motif_wo_impute_nci60/noimptuation_motif_nci60.RDS")
l2 <- do.call(rbind, l2)
l2 <- l2[which(l2$gw.fg.count > 2), ]
l2 <- annotKnownMotifs(l2)
l2$panel <- "NCI60"

ll <- rbind(l2, l1)

ll$moduleColor <- str_split_fixed(rownames(ll), "_|\\.", 3)[, 2]
ll <- ll[, c(18, 17, 1:16)]
write.table(ll, file = "Res/20180314_wgcnaSupTables/moduleMotif.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

