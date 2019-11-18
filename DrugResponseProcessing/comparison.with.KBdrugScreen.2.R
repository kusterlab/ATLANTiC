library(rcellminerData)
library(rcellminer)
library(openxlsx)


kb <- read.xlsx("/media/general/projects/KinobeadDrugScreen/DrugScreen_paper/FINAL supplemental tables for upload/Supplementary Table 1 Inhibitors.xlsx", 
                sheet = 2)

drugMat <- exprs(getRepeatAct(drugData))
annot  <- as(featureData(getRepeatAct(drugData)), "data.frame")
annot.short <- as(featureData(getAct(drugData)), "data.frame")


drug.nop <- tolower(unlist(lapply(annot.short$NAME, function(x) gsub("[[:punct:]]| ", "", x))))
kb.nop <- tolower(unlist(lapply(kb$Drug, function(x) gsub("[[:punct:]]| ", "", x))))

# map by name
intersect(drug.nop, kb.nop)


# map by synonyms
synonyms <- kb$Synonyms
syn <- lapply(strsplit(synonyms, ","), function(x) tolower(gsub("[[:punct:]]| ", "", x)))
tt <- sapply(syn, function(x) x %in% drug.nop)
names(tt) <- names(syn) <- kb$Drug
m2 <- na.omit(sapply(syn[sapply(tt, function(x) any(x))], function(x) x[x %in% c(drug.nop)]))
cbind(rownames(m2), m2)


# map by SMILES
sml <- intersect(annot.short$SMILES, kb$SMILES.code)
s1 <- annot.short[annot.short$SMILES %in% sml, ]
s2 <- kb[kb$SMILES.code %in% sml, ]


fl <- c(intersect(drug.nop, kb.nop), 
        tolower(names(m2)),
        tolower(s2$Drug))
ALL <- unique(fl)

res <- data.frame(
  drug = ALL,
  byName = ALL %in% intersect(drug.nop, kb.nop), 
  bySynonyms = ALL %in% tolower(names(m2)), 
  bySMILES = ALL %in% tolower(s2$Drug), 
  row.names = ALL)

write.table(res, file = "Res/20170329_compare.DTP.KBscreen2/compare2.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

##
kb[kb.nop == "osi027", "SMILES.code"]
annot.short[which(drug.nop == "osi027"), "SMILES"]
