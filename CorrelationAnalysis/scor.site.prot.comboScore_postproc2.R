ll <- readRDS(file = "Res/20170801_comboScoreFittings/summary.RDS")
bm <- do.call(rbind, ll)

any(bm$enzyme == "")
any(is.na(bm$enzyme == ""))

any(bm$entity == "")
any(is.na(bm$entity == ""))

unique(bm$entity)
unique(bm$enzyme)


# top 500

head(bm)
smallestp <- pmin(bm$pval.impute, bm$pval.noimpute, na.rm = TRUE)
fac <- paste(bm$drug1, bm$conc1, bm$drug2, bm$conc2, sep = "__")
resvec <- rep(NA, length(fac))

count <- 0
for (i in unique(fac)) {
  print(count)
  idx <- fac == i
  rk <- rank(smallestp[idx], na.last = TRUE)
  resvec[idx] <- rk <= 500
  count <- count + 1
}

any(is.na(resvec))
table(resvec)

res <- bm[ resvec, ]
dim(res)
head(res)
pm <- pmin(res$pval.impute, res$pval.impute, na.rm = TRUE)
res <- res[!is.na(pm), ]
pm <- pmin(res$pval.impute, res$pval.impute, na.rm = TRUE)
res <- res[pm < 0.01, ]
dim(res)
rownames(res) <- NULL
saveRDS(res, file = "Res/20170801_comboScoreFittings/summary.filtered.RDS")


