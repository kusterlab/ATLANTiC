library(stringr)

files <- c("Res/20170801_comboScoreFittings/scor.imputed.comboScore_site.tryp.RDS", 
           "Res/20170801_comboScoreFittings/scor.imputed.comboScore_site.gluc.RDS", 
           "Res/20170801_comboScoreFittings/scor.imputed.comboScore_prot.tryp.RDS",
           "Res/20170801_comboScoreFittings/scor.imputed.comboScore_prot.gluc.RDS",
           "Res/20170801_comboScoreFittings/scor.noimputation.comboScore_site.tryp.RDS",
           "Res/20170801_comboScoreFittings/scor.noimputation.comboScore_site.gluc.RDS",
           "Res/20170801_comboScoreFittings/scor.noimputation.comboScore_prot.tryp.RDS",
           "Res/20170801_comboScoreFittings/scor.noimputation.comboScore_prot.gluc.RDS")



res <- lapply(files, function(x) {
  print(x)
  xx <- readRDS(file = x)
  xv2 <- str_split_fixed(xx$Var2, "_", n = 4)
  ii <- TRUE
  enz <- ifelse(grepl("tryp", x), "TRYP", "GLUC")
  ent <- ifelse(grepl("site", x), "PHOS", "PROT")
  ds <- data.frame("feature" = xx$Var1[ii], 
                   "drug1" = xv2[ii, 1], 
                   "conc1" = as.numeric(xv2[ii, 2]), 
                   "drug2" = xv2[ii, 3],
                   "conc2" = as.numeric(xv2[ii, 4]), 
                   "cor" = xx$value[ii],
                   "pval" = xx$p[ii],
                   "logp" = -log10(xx$p[ii]),
                   "df" = xx$n[ii]-2,
                   "enzyme"= enz,
                   "entity" = ent
  )
})

library(fastmatch)
library(stringr)

nres2 <- lapply(res, function(x) {
  print("a")
  ids <- paste(x$feature, x$drug1, x$conc1, x$drug2, x$conc2, x$enzyme, x$entity, sep = "_")
  ii <- !grepl("^NA__NA__NA", ids)
  xx <- x[ii, ]
  rownames(xx) <- ids[ii]
  xx
})
#
saveRDS(nres2, file = "Res/20170801_comboScoreFittings/nres2.backup.RDS")

ll <- lapply(1:4, function(i) {
  print(i)
  x1 <- nres2[[i]]
  x2 <- nres2[[i+4]]
  
  cat("correcting rownames \n")
  rownames(x1) <- gsub("935 ERVK-19;ERVK-21;ERVK-6;ERVK-25;ERVK-8;HERVK_113;ERVK-9", 
                       "935 ERVK-19;ERVK-21;ERVK-6;ERVK-25;ERVK-8;HERVK-113;ERVK-9", rownames(x1))
  rownames(x2) <- gsub("935 ERVK-19;ERVK-21;ERVK-6;ERVK-25;ERVK-8;HERVK_113;ERVK-9", 
                       "935 ERVK-19;ERVK-21;ERVK-6;ERVK-25;ERVK-8;HERVK-113;ERVK-9", rownames(x2))
  x1 <- x1[!grepl("^NA.", rownames(x1)), ]
  x2 <- x2[!grepl("^NA.", rownames(x2)), ]
  
  cat("matching union ids \n")
  allid <- unique(c(rownames(x1), rownames(x2)))
  i1 <- fmatch(allid, rownames(x1))
  i2 <- fmatch(allid, rownames(x2))
  
  cat("subsetting matrices \n")
  d1 <- x1[i1, c("cor", "pval", "logp", "df")]
  d2 <- x2[i2, c("cor", "pval", "logp", "df")]
  colnames(d1) <- paste(colnames(d1), "impute", sep = ".")
  colnames(d2) <- paste(colnames(d2), "noimpute", sep = ".")
  
  cat("Binding values \n")
  dc <- cbind(d1, d2)
  
  cat("Splitting meta \n")
  ac <- str_split_fixed(allid, "_", 7)
  colnames(ac) <- c("feature", "drug1", "conc1", "drug2", "conc2", "enzyme", "entity")
  
  cat("binding final \n")
  dcc <- cbind(ac, dc)
  cat("converting to df \n")
  dcc <- data.frame(dcc, stringsAsFactors = FALSE)
  dcc$conc1 <- as.numeric(as.character(dcc$conc1))
  dcc$conc2 <- as.numeric(as.character(dcc$conc2))
  dcc
})


saveRDS(ll, file = "Res/20170801_comboScoreFittings/summary.RDS")

#
bm <- do.call(rbind, ll)
unique(bm$enzyme)
unique(bm$entity)
write.table(bm, file = "Res/20170801_comboScoreFittings/summary.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")






