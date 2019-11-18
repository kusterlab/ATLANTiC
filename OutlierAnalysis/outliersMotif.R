library(PTMotif)

tab <- read.delim(file = "Res/20180328_outlierAnalysis/outlierTable.txt", stringsAsFactors = FALSE)
set <- paste(tab$cell, tab$panel, sep = "_")
slist <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")

#
l1 <- lapply(1:ncol(slist$crc65_expr), function(x) {
  i <- !is.na(slist$crc65_expr[, x])
  gene <- slist$annot$`Gene name`[i]
  gene <- sapply(strsplit(gene, ";"), "[", 1)
  sw <- slist$annot$`Sequence window`[i]
  sw <- sapply(strsplit(sw, ";"), "[", 1)
  df <- unique(data.frame(genes = gene, seqs = sw, stringsAsFactors = FALSE))
  df [nchar(df$seqs) == 31, ]
})
names(l1) <- colnames(slist$crc65_expr)
l2 <- lapply(1:ncol(slist$nci60_expr), function(x) {
  i <- !is.na(slist$nci60_expr[, x])
  gene <- slist$annot$`Gene name`[i]
  gene <- sapply(strsplit(gene, ";"), "[", 1)
  sw <- slist$annot$`Sequence window`[i]
  sw <- sapply(strsplit(sw, ";"), "[", 1)
  unique(data.frame(genes = gene, seqs = sw, stringsAsFactors = FALSE))
  df [nchar(df$seqs) == 31, ]
})
names(l2) <- colnames(slist$nci60_expr)
ll <- c(l1, l2)

all(set %in% names(ll))
length(ll)

lm <- lapply(unique(set), function(s) {
  print(s)
  uu <- function(x) {
    sapply(strsplit(x, ";"), "[", 1)
  }
  i <- set == s
  fs <- uu(tab$Sequence.window[i])
  fg <- uu(tab$Gene.name[i])
  bs <- ll[[s]]$seqs
  bg <- ll[[s]]$genes
  nmin <- max(5, sum(i)*0.05)
  
  dd <- detectMotifs(fg.seqs = fs, bg.seqs = bs, fg.genes = fg, bg.genes = bg,
                     method = "all", ncores = 26, min.seqs = nmin)
  if (is.null(dd))
    return(dd)
  if (nrow(dd) == 0)
    return(dd)
  
  dd$cellline <- s
  dd$overlap <- sapply(dd$motif, function(x) {
    paste(unique(fg[grep(x, fs)]), collapse = ";")
  })
  dd
})
names(lm) <- unique(set)
lm <- do.call(rbind, lm)

write.table(lm, file = "Res/20180328_outlierAnalysis/outlierMotifs.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")










