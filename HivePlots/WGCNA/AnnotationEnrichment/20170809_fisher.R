parseRes <- function(res, subtypes, values, catvec) {
  setnames(res, subtypes)
  res[,values:=values]
  res[,category:=catvec]
  setcolorder(res, c('category', 'values', colnames(res)[!colnames(res)%in%c("category", "values")]))
}

parseQuickGO <- function(res, func_ann, method.all = 'BH', digits = 3, cutoff = 0.05) {
  require(reshape2)
  res <- Reduce(merge, lapply(names(res), function(i) melt(data = res[[i]], id.vars = c('category', 'values'), value.name = i)))
  setkey(res, category)
  res[direction=='more',direction:='enriched']
  res[direction=='less',direction:='depleted']
  res[,`p value`:=signif(`p value`, digits)]
  res[,`q value`:=p.adjust(`q value`, method.all)]
  res <- res[`q value`<=cutoff,]
  if (!is.null(func_ann)) {
    nametable <- func_ann$names
    setkey(nametable, ID)
    res <- nametable[res,,nomatch=NA]
    res[is.na(Name),Name:='']
  }
  return(res)
}

# constructstring <- function(x) {return(paste0(c(paste0('^', x, '$'), paste0('_', x, '$'), paste0('^', x, '_'), paste0('_', x, '_')), collapse = '|'))}
constructstring <- function(x) {return(paste0('(^|_)(', x, ')+(_|$)'))}

fisher <- function(subtypename, category, data, alternative, method.subtype = "BH", digits = 3, multimembersep = NULL, subtypestoexclude = NULL, go_rollup = NULL) {
  require(data.table)
  # require(stringr)
  # For GO-terms, roll up children to current parent category
  if (!is.null(go_rollup)) {
    # cat('Rolled-up GO-terms!\n')
    tmp <- data.frame(as.logical(rowSums(data[,colnames(data)%in%c(category,go_rollup[[category]]),drop=F])))
    colnames(tmp) <- category
    data <- cbind(data[,subtypename,drop=F],tmp)
  }
  len <- length(sort(unique(data[,category])[!is.na(unique(data[,category]))]))
  if (len<=2) {
    values <- tail(sort(unique(data[,category])[!is.na(unique(data[,category]))]), 1)
  } else {
    values <- sort(unique(data[,category])[!is.na(unique(data[,category]))])
  }
  if (!is.null(multimembersep)) {
    subtypes <- as.character(sort(unique(unlist(strsplit(unique(data[,subtypename]), multimembersep))))[!sort(unique(unlist(strsplit(unique(data[,subtypename]), multimembersep))))%in%subtypestoexclude])
    # ymat <- sapply(subtypes, function(i) sapply(values, function(j) sapply(strsplit(data[,subtypename], '_'), function(k) any(k%in%i))))
    res <- lapply(subtypes, function(i) lapply(values, function(j) try(fisher.test(x = data[,category]==j, y = grepl(constructstring(i), data[,subtypename]), alternative = alternative), silent=T)))
    res1 <- lapply(subtypes, function(i) lapply(values, function(j) rownames(data)[grepl(constructstring(i), data[,subtypename])]))
    res2 <- lapply(subtypes, function(i) lapply(values, function(j) rownames(data)[data[,category]==j]))
    res3 <- lapply(subtypes, function(i) lapply(values, function(j) intersect(rownames(data)[data[,category]==j], rownames(data)[grepl(constructstring(i), data[,subtypename])])))
  } else {
    subtypes <- as.character(sort(unique(data[,subtypename]))[!sort(unique(data[,subtypename]))%in%subtypestoexclude])
    res <- lapply(subtypes, function(i) lapply(values, function(j) try(fisher.test(x = data[,category]==j, y = data[,subtypename]==i, alternative = alternative), silent=T))) 
    res1 <- lapply(subtypes, function(i) lapply(values, function(j) rownames(data)[data[,subtypename]==i]))
    res2 <- lapply(subtypes, function(i) lapply(values, function(j) rownames(data)[data[,category]==j]))
    res3 <- lapply(subtypes, function(i) lapply(values, function(j) intersect(rownames(data)[data[,category]==j], rownames(data)[data[,subtypename]==i])))
  }
  nareplace <- lapply(res, function(i) which(unlist(lapply(i, function(j) class(j)=="try-error"))))
  whichreplace <- which(lapply(nareplace, length)>0)
  for (i in whichreplace) {
    tmp <- fisher.test(x = rep(c(T,T,F), times=4), y=rep(c(T,F,F), times=4), alternative = alternative)
    tmp$p.value <- NA
    tmp$conf.int <- NA
    tmp$estimate <- NA
    tmp$null.value <- NA
    res[[i]][nareplace[[i]]] <- lapply(seq_along(nareplace[[i]]), function(j) tmp)
  }
  catvec <- rep(category, times=length(values))
  sublength <- data.table(do.call(cbind, lapply(res1, function(i) try(sapply(i, function(j) length(j)), silent=T))))
  term <- data.table(do.call(cbind, lapply(res2, function(i) try(sapply(i, function(j) paste0(j, collapse = ';')), silent=T))))
  termlength <- data.table(do.call(cbind, lapply(res2, function(i) try(sapply(i, function(j) length(j)), silent=T))))
  inter <- data.table(do.call(cbind, lapply(res3, function(i) try(sapply(i, function(j) paste0(j, collapse = ';')), silent=T))))
  interlength <- data.table(do.call(cbind, lapply(res3, function(i) try(sapply(i, function(j) length(j)), silent=T))))
  univlength <- data.table(matrix(rep(nrow(data), times = length(subtypes)*length(values)), ncol = length(subtypes), byrow = F))
  expeclength <- data.table(matrix(unlist(termlength)/unlist(univlength)*unlist(sublength), ncol = length(subtypes), byrow = F))
  pvals <- data.table(do.call(cbind, lapply(res, function(i) try(sapply(i, function(j) j$p.value), silent=T))))
  or <- data.table(signif(do.call(cbind, lapply(res, function(i) sapply(i, function(j) j$estimate))), digits = digits))
  dir <- data.table(ifelse(or>1, "more", "less"))
  
  sublength <- parseRes(sublength, subtypes, values, catvec)
  term <- parseRes(term, subtypes, values, catvec)
  termlength <- parseRes(termlength, subtypes, values, catvec)
  inter <- parseRes(inter, subtypes, values, catvec)
  interlength <- parseRes(interlength, subtypes, values, catvec)
  univlength <- parseRes(univlength, subtypes, values, catvec)
  expeclength <- parseRes(expeclength, subtypes, values, catvec)
  pvals <- parseRes(pvals, subtypes, values, catvec)
  or <- parseRes(or, subtypes, values, catvec)
  dir <- parseRes(dir, subtypes, values, catvec)
  
  result <- list(p.values=pvals,
                 `odds ratio`=or,
                 `direction`=dir,
                 `count in subtype`=sublength,
                 `in term`=term,
                 `count in term`=termlength,
                 `overlap`=inter,
                 `count in overlap`=interlength,
                 `count in universe`=univlength,
                 `expected`=expeclength)
  return(result)
}

padj <- function(pvals, method.subtype = "BH", digits = 3) {
  fpvals <- as.matrix(pvals[,!colnames(pvals)%in%c("category", "values"),with=F])
  res <- matrix(apply(fpvals, 1, function(i) signif(p.adjust(i, method=method.subtype, n = length(i)), digits = digits)), nrow = nrow(fpvals), ncol = ncol(fpvals), byrow = T)
  # res <- signif(matrix(p.adjust(c(fpvals), method=method), nrow=nrow(fpvals), ncol=ncol(fpvals),byrow=F), digits = digits)
  colnames(res) <- colnames(pvals)[!colnames(pvals)%in%c("category", "values")]
  res <- data.table(pvals[,colnames(pvals)%in%c("category", "values"),with=F], res, check.names = F)
  return(res)
}

fisher.rbind <- function(subtypename, categories, data, parseResult = F, func_ann = NULL, alternative="two.sided", method.subtype = "none", method.all = 'BH', digits = 3, cutoff = 0.05, multimembersep = NULL, subtypestoexclude = NULL, go_rollup = NULL, ncores = 10) {
  require(data.table)
  require(pbmcapply)
  if (!is.null(multimembersep)) {
    cat(paste0("Subtypes are split at '", multimembersep, "'\n"))
  }
  data[data=="NA"] <- NA
  # data <- data[,c(subtypename, categories)]
  categories <- names(which(colSums(is.na(data))!=nrow(data)))[names(which(colSums(is.na(data))!=nrow(data)))%in%categories]
  tmp <- data[,colSums(is.na(data))!=nrow(data)]
  if(ncol(tmp)<ncol(data)) {
    cat("Some columns only contain NAs - dropped!")
  }
  data <- tmp
  res <- pbmclapply(colnames(data)[match(categories, colnames(data))], function(i) {
      fisher(subtypename = subtypename, category = i, subtypestoexclude = subtypestoexclude, data = data, alternative = alternative, multimembersep = multimembersep, go_rollup = go_rollup)
    }, mc.cores = ncores)
  names(res) <- colnames(data)[match(categories, colnames(data))]
  pvals <- rbindlist(lapply(res, function(i) i$p.values))
  or <- rbindlist(lapply(res, function(i) i$`odds ratio`))
  dir <- rbindlist(lapply(res, function(i) i$direction))
  qvals <- padj(pvals, method = method.subtype, digits = digits)
  sublength <- rbindlist(lapply(res, function(i) i$`count in subtype`))
  term <- rbindlist(lapply(res, function(i) i$`in term`))
  termlength <- rbindlist(lapply(res, function(i) i$`count in term`))
  inter <- rbindlist(lapply(res, function(i) i$`overlap`))
  interlength <- rbindlist(lapply(res, function(i) i$`count in overlap`))
  univlength <- rbindlist(lapply(res, function(i) i$`count in universe`))
  expeclength <- rbindlist(lapply(res, function(i) i$`expected`))
  result <- list(
    `count in universe`=univlength,
    `count in term`=termlength,
    `count in variable`=sublength,
    `count in overlap`=interlength,
    `expected count`=padj(pvals = expeclength, method.subtype = 'none', digits = digits),
    `direction`=dir,
    `p value`=pvals,
    `q value`=qvals,
    `odds ratio`=or,
    `in term`=term,
    `overlap`=inter
    )
  if(parseResult) {
    result <- parseQuickGO(result, func_ann, method.all, digits, cutoff)
  }
  return(result)
}