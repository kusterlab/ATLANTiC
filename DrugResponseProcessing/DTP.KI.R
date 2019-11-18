library(XML)
library(openxlsx)

dtp2 <- readHTMLTable("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Dat/DTP_20170328/DTP_NCI60_ZSCORE.html", header = TRUE)
kb <- read.xlsx("/media/general/projects/NCI_60_phospho/Chen/Dat/kbscreen&review/Supplementary Table 1 Inhibitors_withcitations.xlsx", sheet = 2)
bkc <- read.delim("/media/general/projects/NCI_60_phospho/Chen/Dat/kbscreen&review/kinobeadScreen.txt")
idmap <- read.delim("/media/general/projects/NCI_60_phospho/Chen/Dat/kbscreen&review/bkCID2SID.txt")
expr <- readRDS("Res/20170719_nci60.proteinGroups.gluc/nci60.proteinGroups.gluc.RDS")

# s1
s1 <- dtp2[grepl("K", dtp2$`Mechanism of action`)|grepl("nib$", dtp2$`Drug name`), c("Drug name", "Mechanism of action")]

# s3
s3 <- dtp2[dtp2$`PubChem SID` %in% idmap$SID, c("Drug name", "Mechanism of action")]
setdiff(s3$`Drug name`, s1$`Drug name`)

# s4
s4 <- na.omit(dtp2[dtp2$SMILES %in% kb$SMILES.code, c("Drug name", "Mechanism of action")])
setdiff(s4$`Drug name`, s1$`Drug name`)
setdiff(s3$`Drug name`, s4$`Drug name`)


# 
colnames(dtp2)
sel <- dtp2[dtp2$`Drug name` %in% unique(c(as.character(s1$`Drug name`), as.character(s3$`Drug name`), as.character(s4$`Drug name`))), 
            c(1:6, 67:68)]
i <- match(sel$SMILES, kb$SMILES.code)
sel <- cbind(sel, kb[i, c(1, 5:8)])

sel <- within(sel, {
  dn <- as.character(get("Drug name"))
  udn <- unique(dn)
  expa <- function(x) {
    v <- na.omit(x)
    if (length(v) > 0)
      x[is.na(x)] <- v
    return(x)
  }
  for (i in udn) {
    ii <- dn == i
    Drug[ii] <- expa(Drug[ii])
    Designated.targets[ii] <- expa(Designated.targets[ii])
    Other.targets[ii] <- expa(Other.targets[ii])
    Binding.mode[ii] <- expa(Binding.mode[ii])
    Binding.type[ii] <- expa(Binding.type[ii])
  }
  rm(list = c("expa", "dn", "udn", "i", "ii"))
})



f00 <- dtp2[rownames(sel), ]
identical(f00$`Drug name`, sel$`Drug name`)

dtp.gi50 <- apply(dtp2[rownames(sel), -c(1:6, 67:68)], 2, as.numeric)


te <- gsub("-|/|\\(TB\\)|/ATCC| ", "_", colnames(dtp.gi50))
te <- gsub("T_47D", "T47D", te)
te <- gsub("HS_578T", "HS578T", te)
te <- gsub("LOX_IMVI", "LOXIMVI", te)
te <- gsub("A549_", "A549", te)
te <- gsub("HL_60_", "HL_60", te)
te <- gsub("COLO_205", "COLO205", te)


colnames(dtp.gi50) <- te

setdiff(colnames(dtp.gi50), expr$xref$name_DTP)
setdiff(expr$xref$name_DTP, colnames(dtp.gi50))


dtp.gi50 <- dtp.gi50[, match(expr$xre$name_DTP, colnames(dtp.gi50))]

cbind(colnames(dtp.gi50), expr$xref$name_MQ)
colnames(dtp.gi50) <- expr$xref$name_AUC


dim(dtp.gi50)
dim(sel)


res <- list(annot = sel, gi50 = dtp.gi50)
saveRDS(res, file = "Res/20170907_DTPKIprocess/dtpki.RDS")
