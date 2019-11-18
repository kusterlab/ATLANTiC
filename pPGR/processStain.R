library(openxlsx)
library(gdata)
library(readxl)
library(reshape2)
library(beeswarm)

ll <- lapply(1:17, function(s) read_excel(path = "Dat/PGR/combined.xlsx", sheet = s))

#
array <- lapply(1:34, function(s) as.matrix(read_excel(path = "Dat/PGR/array.xlsx", sheet = s)))
array <- sapply(array, function(x) {
  if (sum(dim(x)) == 0)
    return(NA)
  melt(x)
  })

parsestain <- function(idx, stain, i) {
  
  stain <- stain[!is.na(stain$`H-Nr. , Block`) & stain$`H-Nr. , Block` != "Pla", ]
  perc <- cbind(array[[idx[1]]]["value"], array[[idx[2]]]["value"])
  colnames(perc) <- c("id", "value")
  perc$id <- as.character(perc$id)
  perc$value <- as.character(perc$value)
  perc <- perc[ perc$id %in% stain$`H-Nr. , Block`, ]
  perc <- perc[order(perc$id), ]
  stain <- stain[order(stain$`H-Nr. , Block`), ]
  
  perc$value[perc$value == "n.a."] <- NA
  perc$value[perc$value == "kT"] <- 0
  perc$value[perc$value == "+"] <- 0.6
  perc$value <- as.numeric(perc$value)
  
  if (!all(perc$id == stain$`H-Nr. , Block`))
    stop("patient id different")
  
  
  intens <- grep("intens", colnames(stain))
  percCol <- grep("%", colnames(stain))
  res <- data.frame(
    # Int.Nr = stain$`Int. Nr.`, 
    HNR.block = stain$`H-Nr. , Block`, 
    percent.tumor = perc$value,
    stain.percent.protein = as.numeric(stain[[percCol[1]]]),
    stain.intensity.protein = stain[[intens[1]]],
    stain.percent.phos = as.numeric(stain[[percCol[2]]]),
    stain.intensity.phos = stain[[intens[2]]],
    silver = i
  )
}


n <- 1:17
mat <- matrix(1:34, ncol = 2, byrow = TRUE)
stainList <- lapply(n, function(i) {
  print(i)
  parsestain(idx = mat[i, ], stain = ll[[i]], i = i)
  })
stained <- do.call(rbind, stainList)



saveRDS(stained, file = "Res/20170707_pGR/pgrstain.RDS")
# stained <- unique(stained)


##
stained <- readRDS(file = "Res/20170707_pGR/pgrstain.RDS")


png("Res/20170707_pGR/pgrstain_percent_intensity.png", width=32, height = 12, unit = "cm", res = 300)
layout(matrix(1:2, 1, 2))
m1 <- table(stained$stain.percent.protein, stained$stain.intensity.protein)
barplot(t(m1), beside = TRUE, legend.text = colnames(m1), xlab = "% of staining cells", ylab = "count", main = "protein",
        col = rev(gray.colors(4)), args.legend = list(x = "topright", bty = "n", title = "Intensity"))

m1 <- table(stained$stain.percent.phos, stained$stain.intensity.phos)
barplot(t(m1), beside = TRUE, legend.text = colnames(m1), xlab = "% of staining cells", ylab = "count", main = "phospho",
        col = rev(gray.colors(5)), args.legend = list(x = "topleft", bty = "n", title = "Intensity", title.adj = 0.1))
dev.off()

