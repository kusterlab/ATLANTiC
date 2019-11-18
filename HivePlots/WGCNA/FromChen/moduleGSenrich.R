library(PTMotif)
library(stringr)
library(data.table)


crc <- fread("/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/res/20180302_wgcna_allsites_allproteins_incnas_crc65_nonas_crc65/CSVs/GO/goall_NoImputation.csv", 
             integer64 = "double", data.table = FALSE)

nci <- fread("/media/kusterlab/users_files/martin/postdoc/projects/phosphoproject/res/20180308_wgcna_allsites_allproteins_incnas_nci60_nonas_nci60/CSVs/GO/goall_NoImputation_nci60.csv", 
             integer64 = "double", data.table = FALSE)


crc$panel <- "CRC65"
nci$panel <- "NCI60"

df <- rbind(crc, nci)
df <- df[c(4, 1:3, 5:17)]
head(df)

df <- df[order(df$`q value`, decreasing = FALSE), ]
df <- df[order(df$variable), ]
df <- df[order(df$panel), ]

write.table(df, file = "Res/20180314_wgcnaSupTables/moduleGSenrich.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

