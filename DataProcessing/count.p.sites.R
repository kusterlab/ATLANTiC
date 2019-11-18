library(data.table)
library(matrixStats)
source("../R/impute2.R")

sites <- fread("../Final_dataset/MQ1.5.5.1_pNCI60_pCRC65_complete/txt/Phospho (STY)Sites.txt", 
               integer64 = "double", data.table = FALSE)
table(sapply(sites, class))

# =========================== expression matrix ================================
ii <- grepl("Intensity ", colnames(sites)) & !grepl("___", colnames(sites))
table(ii)

expr <- apply(sites[, ii], 2, as.numeric)
colnames(expr) <- gsub("Intensity ", "", colnames(expr))

# two panels
i_crc <- grepl("CRC", colnames(expr))
i_nci <- grepl("NCI60", colnames(expr))

# filtering
f_rev <- !grepl("^REV_", sites$`Leading proteins`)
f_only0 <- rowSums(expr) > 0
f_loc <- sites$`Localization prob` >= 0.75
f_nomsms <- sites$`MS/MS IDs` != ""
f_all <- f_rev & f_only0 & f_loc & f_nomsms

# transform and filtering
expr <- log10(expr)
expr[is.infinite(expr)] <- NA
expr <- sweep(expr, 2, colMedians(expr, na.rm = TRUE), "-") + median(expr, na.rm = TRUE)
expr <- expr[f_all, ]


ss <- paste(sapply(strsplit(sites$Protein, "-|;"), "[", 1), sites$Position, sep = "_")
ss <- ss[f_all]
length(unique(ss))
table(f_all)


###########
sites.gluc <- fread("../Final_dataset/MQ1.5.5.1_pNCI60_GluC_complete/txt/Phospho (STY)Sites.txt", 
               integer64 = "double", data.table = FALSE)
table(sapply(sites.gluc, class))

# =========================== expression matrix ================================
ii <- grepl("Intensity ", colnames(sites.gluc)) & !grepl("___", colnames(sites.gluc))
table(ii)

expr <- apply(sites.gluc[, ii], 2, as.numeric)
colnames(expr) <- gsub("Intensity ", "", colnames(expr))

# two panels
i_crc <- grepl("CRC", colnames(expr))
i_nci <- grepl("NCI60", colnames(expr))

# filtering
f_rev <- !grepl("^REV_", sites.gluc$`Leading proteins`)
f_only0 <- rowSums(expr) > 0
f_loc <- sites.gluc$`Localization prob` >= 0.75
f_nomsms <- sites.gluc$`MS/MS IDs` != ""
f_all <- f_rev & f_only0 & f_loc & f_nomsms

# transform and filtering
expr <- log10(expr)
expr[is.infinite(expr)] <- NA
expr <- sweep(expr, 2, colMedians(expr, na.rm = TRUE), "-") + median(expr, na.rm = TRUE)
expr <- expr[f_all, ]

ss2 <- paste(sapply(strsplit(sites.gluc$Protein, "-|;"), "[", 1), sites.gluc$Position, sep = "_")
ss2 <- ss2[f_all]
length(ss2)

length(unique(ss))
length(unique(ss2))
length(unique(c(ss, ss2)))



