library(RTCGAToolbox)
library(genefilter)
library(TCGAbiolinks)

pat <- readRDS(file = "Res/20170324_Fluororacil_TCGA/tcga.patient.drug.RDS")

bc <- split(pat$bcr_patient_barcode, pat$CancerType)

getFirehoseAnalyzeDates()
names(bc)

mrnas <- lapply(names(bc), function(tc) {
  print(tc)
  t1 <- getFirehoseData(
    tc, runDate = "20160128", gistic2_Date = NULL,
    RNAseq_Gene = FALSE, Clinic = FALSE, miRNASeq_Gene = FALSE,
    RNAseq2_Gene_Norm = TRUE, CNA_SNP = FALSE, CNV_SNP = FALSE,
    CNA_Seq = FALSE, CNA_CGH = FALSE, Methylation = FALSE,
    Mutation = FALSE, mRNA_Array = FALSE, miRNA_Array = FALSE,
    RPPA = FALSE, RNAseqNorm = "raw_counts",
    RNAseq2Norm = "normalized_count", 
    forceDownload = FALSE, destdir = ".",
    fileSizeLimit = 500, 
    getUUIDs = FALSE)
  t1@RNASeq2GeneNorm
})
names(mrnas) <- names(bc)

bm <- do.call(cbind, mrnas)
dim(bm)

colnames(bm) <- substr(colnames(bm), 1, 12)

isub <- intersect(pat$bcr_patient_barcode, colnames(bm))
pat.sub <- pat[pat$bcr_patient_barcode %in% isub, ]
dim(pat.sub)

smat <- bm[, pat.sub$bcr_patient_barcode]
smat <- log10(smat+1)
hist(smat, breaks = 100)


table(pat.sub$measure_of_response)

##
icr <- pat.sub$measure_of_response == "Complete Response"
icpd <- pat.sub$measure_of_response == "Clinical Progressive Disease"
imid <- pat.sub$measure_of_response == "Partial Response" | pat.sub$measure_of_response == "Stable Disease"


rt <- rowttests(cbind(smat[, icr], smat[, imid|icpd]), fac = factor(rep(c("CR", "NCR"), c(sum(icr), sum(imid|icpd)))))
rt <- rt[order(rt$p.value, decreasing = FALSE), ]
rt$fdr <- p.adjust(rt$p.value, method = "fdr")


beeswarm(smat["SLC9A3", ] ~ icr)
boxplot(smat["SLC9A3", ] ~ icr)
t.test(smat["SLC9A3", ] ~ icr, var.equal = TRUE)





# compare icr and icpd
rt1 <- rowttests(cbind(smat[, icr], smat[, imid]), fac = factor(rep(c("CR", "MR"), c(sum(icr), sum(imid)))))
rt2 <- rowttests(cbind(smat[, imid], smat[, icpd]), fac = factor(rep(c("MR", "CPD"), c(sum(imid), sum(icpd)))))

if1 <- rt1$dm > 0 & rt1$p.value < 0.05  & rt2$dm < 0 & rt2$p.value < 0.05
if2 <- rt1$dm < 0 & rt1$p.value < 0.05  & rt2$dm > 0 & rt2$p.value < 0.05

mat1 <- cbind(rt1, rt2)[if1, ] # sensitive markers
mat1 <- mat1[order(mat1[, 6] * mat1[, 3], decreasing = FALSE), ]

mat2 <- cbind(rt1, rt2)[if2, ] # resistance marker
mat2 <- mat2[order(mat2[, 6] * mat2[, 3], decreasing = FALSE), ]

write.table(mat1, file = "Res/20170324_Fluororacil_TCGA/sensitive.marker.tcga.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(mat2, file = "Res/20170324_Fluororacil_TCGA/resistance.marker.tcga.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

################################################################################

rt3 <- rowttests(cbind(smat[, icr], smat[, icpd]), fac = factor(rep(c("CR", "CPD"), c(sum(icr), sum(icpd)))))

markers <- c("ABCC3", "ABCC4", "TK1", "TYMS", "UMPS", "TYMP", "RRM1", "RRM2", "PPAT", 
             "DPYD", "UPP1", "UCK2", "UCK1", "SLC29A1")
rt1[markers, ]
rt2[markers, ]
rt3[markers, ]

library(beeswarm)
beeswarm(smat["RAD52", ] ~ as.character(pat.sub$measure_of_response))
boxplot(smat["RAD52", ] ~ pat.sub$measure_of_response)

identical(names(smat["RAD52", ]), 
      pat.sub$bcr_patient_barcode)

t.test(smat["RAD52", icr], smat["RAD52", icpd])

t.test(smat["RAD52", icr], smat["RAD52", imid], var.equal = TRUE)
t.test(smat["RAD52", imid], smat["RAD52", icpd], var.equal = TRUE)
