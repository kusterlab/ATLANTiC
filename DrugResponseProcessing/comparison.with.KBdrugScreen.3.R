
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
si <- intersect(drug.nop, kb.nop)
tab1 <- annot.short[drug.nop %in% si, ]


# map by SMILES
sml <- intersect(annot.short$SMILES, kb$SMILES.code)
tab2 <- annot.short[annot.short$SMILES %in% sml, ]

taba <- unique(rbind(tab1, tab2))


ss <- annot$nsc %in% taba$NSC
smat <- drugMat[ss, ]
amat <- annot[ss, ]
vec <- structure(annot.short$NAME, names = annot.short$NSC)
amat$drug <- vec[amat$nsc]

dat <- list(annotation = amat, 
            drugAct = smat)
sapply(dat, dim)


saveRDS(dat, "Res/20170329_compare.DTP.KBscreen2/measuredKIs.dtp.RDS")

png("Res/20170329_compare.DTP.KBscreen2/conc.KI.dtp.png", width = 28, height = 12, units = "cm", res = 150)
par(mar = c(8, 5, 0.5, 0.5))
ord <- order(amat$drug)
boxplot(t(smat)[, ord], las = 2, names = amat$drug[ord], col = c(amat$used_in_zscore)+2,
        ylab = "-log10(GI50Molar)", cex = 0.5, cex.axis = 0.6)
legend("bottomleft", col = 2:3, pch = 15, bty = "n", 
       legend = c("QC not passed", "QC passed"))
legend("bottomright", col = c("orange", "cyan"), lty=1, bty = "n", 
       legend = c("3 microM", "10 microM"))

# text(x=2.5, y=5, "3 microM")
abline(h = -log10(3*10^-6), col = "orange")
abline(h = -log10(10*10^-6), col = "cyan")
# text(x=2.5, y=4, "3 microM")
dev.off()

head(amat)
dim(amat)
annot.short[ which(annot.short$NAME == "Fluorouracil"), ]


