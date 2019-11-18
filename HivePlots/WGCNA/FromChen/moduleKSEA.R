library(PTMotif)
library(stringr)

ksdb <- readRDS("Res/20160803_kinaseSubstrateDatabaseParse/integDatabase.RDS")
ksdb$label <- paste(ksdb$sub_geneName, ksdb$sub_seq7)
ksdb$set <- paste(ksdb$kinase, ksdb$database)

mod <- read.delim(file = "Res/20180314_wgcnaSupTables/moduleCandidates.txt", header = TRUE, stringsAsFactors = FALSE)
mod <- mod[!is.na(mod$sequenceWindow), ]
mod$label <- paste(mod$gene, substr(mod$sequenceWindow, 9, 23))
mod$set <- paste(mod$moduleColor, mod$panel)

iit <- intersect(mod$label, ksdb$label)

ksdb <- ksdb[ksdb$label %in% iit, ]
ksdb <- ksdb[ksdb$set %in% names(table(ksdb$set)[table(ksdb$set) > 9]), ]
ksdb$set <- as.factor(ksdb$set)

mod <- mod[mod$label %in% iit, ]
mod <- mod[mod$set %in% names(table(mod$set)[table(mod$set) > 9]), ]

##
setn <- table(ksdb$set) # white ball
totn <- length(unique(mod$label))
setno <- totn - setn # black ball
qq <- tapply(paste(ksdb$sub_geneName, ksdb$aminoAcid, ksdb$position, ksdb$sub_seq7), ksdb$set, 
       paste, collapse = ";")


modl <- split(mod, mod$set)


lm <- lapply(names(modl), function(x) {
  
  x0 <- strsplit(x, ' ')[[1]]
  moduleColor <- x0[1]
  pnl<- x0[2]
  m0 <- modl[[x]]
  kselect <- nrow(m0)
  ii <- ksdb$label %in% m0$label
  kssub <- ksdb[ii, ]
  q0 <- tapply(paste(kssub$sub_geneName, kssub$aminoAcid, kssub$position, kssub$sub_seq7), kssub$set, 
               paste, collapse = ";")
  
  nw <- table(ksdb$set[ii])
  
  nn <- names(nw)
  nn <- str_split_fixed(nn, pattern = " ", 2)
  
  n.overlap <- nw
  n.modulenotset <- kselect - n.overlap 
  n.setnotmodule <- setn - n.overlap
  n.notsetmodule <- setno - kselect + n.overlap
  or <- (n.overlap*n.notsetmodule)/(n.modulenotset*n.setnotmodule)
  expect <- round(kselect * (setn/totn), digits = 2)
  
  pv <- phyper(pmax(nw-1, 0), setn, setno, kselect, lower.tail = FALSE)
  fdr <- p.adjust(pv, method = "fdr")
  
  df <- data.frame(
    moduleColor = x0[1],
    kinase = c(nn[, 1]),
    count.universe = c(n.notsetmodule),
    count.targets = c(n.setnotmodule),
    count.module = c(n.modulenotset),
    count.overlap = c(n.overlap),
    expected.count = c(expect),
    p.value = c(pv),
    q.value = c(fdr),
    odds.ratio = c(or),
    targets = c(qq[names(nw)]),
    overlap = c(q0[names(nw)]),
    database = c(nn[, 2]),
    panel = c(pnl),
    stringsAsFactors = FALSE
  )
  df <- df[df$q.value < 0.05, ]
  
})

llm <- do.call(rbind, lm)


write.table(llm, file = "Res/20180314_wgcnaSupTables/moduleUpstreatKinase.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


