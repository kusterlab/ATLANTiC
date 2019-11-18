source("R/kinase_substrate_databaseIntegration/uniprot.XREF.R")
library(seqinr)

# processing of three databases, including: phosphositeplus, phosphoELM and KEA2.0

## PhosphoSitePlus
psp <- read.delim("Dat/phosphoSitePlus/New/Kinase_Substrate_Dataset", comment.char = "#", stringsAsFactors = FALSE)
psp <- psp[psp$KIN_ORGANISM == "human" & psp$SUB_ORGANISM == "human", ]

psp$aminoAcid <- substr(psp$SUB_MOD_RSD, 1, 1)
psp$position <- as.integer(substr(psp$SUB_MOD_RSD, 2, nchar(psp$SUB_MOD_RSD)))
psp$id <- accTo(sapply(strsplit(psp$SUB_ACC_ID, "-"), "[", 1), "UniProtKB-ID")
head(psp)

write.table(psp, file = "Res/20160803_kinaseSubstrateDatabaseParse/phosphoSitePlus.selftab.txt",
            col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")

## Phospho.ELM
pelm <- read.delim("Dat/phosphoELM/phosphoELM_all_2015-04.dump", stringsAsFactors = FALSE)
pelm <- pelm[pelm$species == "Homo sapiens" & pelm$kinases != "", ]

fix <- paste(rep("_", 15), collapse = "")
pelm$sequence <- paste(fix, pelm$sequence, fix, sep="")
seqs <- substr(pelm$sequence, pelm$position, pelm$position+30)

pelm$sequence <- seqs
head(pelm)

pelm$geneName <- accTo(sapply(strsplit(pelm$acc, "-"), "[", 1), "Gene_Name")
pelm$ID <- accTo(sapply(strsplit(pelm$acc, "-"), "[", 1), "UniProtKB-ID")
write.table(pelm, file = "Res/20160803_kinaseSubstrateDatabaseParse/phosphoELM.selftab.txt",
            col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")

## KEA2.0
kea <- readLines("Dat/KEA2.0/kinase-substrate_phospho-site_level_set_library.tsv")
kea <- strsplit(kea, split = "\t")
kinase <- sapply(kea, "[", 1)
kea <- lapply(kea, "[", -(1:2))
names(kea) <- kinase

fa <- read.fasta("Dat/Uniprot/uniprot_sprot.fasta", seqtype = "AA")
x <- unlist(kea)
m <- do.call("rbind", strsplit(x, split = "_"))
aa <- substr(m[, 2], 1, 1)
pos <- as.integer(gsub("-", "", substr(m[, 2], 2, nchar(m[, 2]))))
acc <- gnToAcc(m[, 1])
acc1 <- sapply(strsplit(acc, ";"), "[", 1)
ID <- accTo(acc1, "UniProtKB-ID")
mm <- cbind(rep(kinase, sapply(kea, length)), x, aa, pos, acc1, ID)
mm <- unique(mm)

allnames <- paste("sp", mm[, "acc1"], mm[, "ID"], sep = "|")
allnames[!allnames %in% names(fa)]

ffa <- fa[allnames]
ffa <- sapply(ffa, function(x) {
  paste(c(rep("_", 15), x, rep("_", 15)), collapse = "")
})

flankseq15 <- substr(ffa, as.integer(mm[, "pos"]), as.integer(mm[, "pos"])+30)
aminoAcid <- substr(flankseq15, 16, 16)
barplot(table(aminoAcid))

mm <- cbind(mm, flankseq15, aminoAcid)
head(mm)
table(mm[, "aa"] == mm[, "aminoAcid"])

mm <- mm[mm[, "aa"] == mm[, "aminoAcid"] & nchar(mm[, "flankseq15"]) == 31, ]
write.table(mm, file = "Res/20160803_kinaseSubstrateDatabaseParse/KEA.selftab.txt",
            col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")


