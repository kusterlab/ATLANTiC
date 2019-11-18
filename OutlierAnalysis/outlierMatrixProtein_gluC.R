source("https://raw.githubusercontent.com/mengchen18/RFunctionCollection/master/outliers.R")

library(openxlsx)
library(matrixStats)
if (!dir.exists("Res/20180328_outlierAnalysis"))
  dir.create("Res/20180328_outlierAnalysis")


t1 <- read.xlsx("../Manuscript/current/supplementTables/table_S3_Intensity_of_proteins.xlsx", sheet = 7)
an <- read.xlsx("../Manuscript/current/supplementTables/table_S3_Intensity_of_proteins.xlsx", sheet = 6)
nam <- t1$ID
t1$ID <- NULL
t1 <- apply(t1, 2, as.numeric)
rownames(t1) <- nam

i <- findOutlier(t1, foldthresh = 5, pvalue = 0.1, reachLowBound = FALSE, window = 0.5)

getab <- function(m, a, col = make.names(c("Gene.names", "iBAQ", "Majority protein IDs", "Protein names"))) {
  i <- findOutlier(m, foldthresh = 5, pvalue = 0.1, reachLowBound = FALSE, window = 0.5)
  ls <- lapply(names(i$outlierIndexColumns), function(ii) {
    ir <- i$outlierIndexColumns[[ii]]
    d <- a[ir, col]
    if (nrow(d) == 0)
      return(d)
    d$cellline <- ii
    d
  })
  do.call(rbind, ls)
}

tab <- getab(t1, a = an)
tab$cellline <- paste0(tab$cellline, "_NCI60_GluC")
write.table(tab, file = "Res/20180328_outlierAnalysis/outlierTableProtein_GluC.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
