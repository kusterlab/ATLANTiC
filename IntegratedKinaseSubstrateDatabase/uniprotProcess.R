
up <- readLines("Dat/Uniprot/uniprot_sprot_human.dat")
spl <- cumsum(up == "//")
uniqueEntry <- split(up, spl)

parsePPsites <- function(entry) {
  i1 <- grepl("FT   MOD_RES", entry) & grepl("Phospho", entry)
  i2 <- grepl("by", entry)
  i3 <- grepl("\\.", entry)
  s1 <- which(i1 & i2)
  
  if (length(s1) == 0)
    return(NULL)
  
  s2 <- which(i3)
  ali <- intersect(s1, s2)
  umap <- s1[!s1 %in% s2]
  if (length(umap) > 0) {
    cr <- lapply(umap, function(x) {
      uStart <- x
      uEnd <- s2[match(TRUE, s2 > uStart)]
      entry[uStart:uEnd]
    })
    cr <- sapply(cr, function(x) {
      x[-1] <- gsub("and", ",", trimws(gsub("FT", "", x[-1])))
      paste(x, collapse = "")
    })
  } else
    cr <- NULL
  str <- c(entry[ali], cr)
  str <- gsub(" and| or", ",", str)
  
  
  str <- strsplit(str, ";")
  str <- sapply(str, function(x) {
    x1 <- strsplit(x[1], " ")
    x1 <- x1[[1]][nchar(x1[[1]]) > 0]
    
    aa <- switch(x1[5], 
                 "Phosphotyrosine" = "Y",
                 "Phosphothreonine" = "T",
                 "Phosphoserine" = "S")
    
    x2 <- gsub("by|and|or|\\{(.*)|\\.", ",", x[[2]])
    x2 <- gsub("alternate|in vivo|in vitro", "", x2)
    x2 <- trimws(strsplit(x2, ",")[[1]])
    x2 <- x2[nchar(x2)>0]
    c(x1[3], aa, paste(x2, collapse=";"))
  })
  
  t(str)
}



accs <- sapply(uniqueEntry, function(entry) {
  entry <- entry[entry != "//"]
  if (length(entry) < 2) 
    return(c(NA, NA))
  id <- strsplit(entry[1], " ")
  id <- id[[1]][nchar(id[[1]]) > 0]
  acc <- strsplit(entry[2], "( |;)")
  acc <- acc[[1]][nchar(acc[[1]]) > 0]
  c(id[2], acc[2])
})

sites <- lapply(uniqueEntry, function(entry) {
  entry <- entry[entry != "//"]
  if (length(entry) < 2) 
    return(NULL)
  parsePPsites(entry = entry)
})

names(sites) <- accs[1, ]
sites <- sites[!sapply(sites, is.null)]

smat <- do.call("rbind", lapply(sites, as.data.frame))
colnames(smat) <- c("site", "aminoAcid", "kinase")
smat$protein <- sapply(strsplit(rownames(smat), "\\."), "[", 1)

colnames(accs) <- accs[1, ]
smat$proteinAcc <- accs[2, smat$protein]
smat <- as.matrix(smat)

smatExt <- lapply(1:nrow(smat), function(i) {
  x <- smat[i, ]
  if (grepl(";", x[3])) {
    ki <- strsplit(x[3], ";")$kinase
    m <- matrix(rep(x, length(ki)), nrow = length(ki), byrow = TRUE)
    m[, 3] <- ki
    m
  } else
    m <- x
  m
})

smat2 <- do.call(rbind, smatExt)
colnames(smat2)
write.table(smat2, file = "Res/20160803_kinaseSubstrateDatabaseParse/uniprot.parse.intermid.txt",
            col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)

