# https://stackoverflow.com/questions/36220585/efficient-jaccard-similarity-documenttermmatrix?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
library(data.table)
library(text2vec)
library(magrittr)
library(Matrix)
require(reshape2)

generateTermAssignment <- function(ids, func_ann, identifier = 'SYMBOL') {
  require(pbmcapply)
  setkeyv(func_ann$annotations, identifier)
  nongo <- func_ann$annotations[.(ids),,nomatch = 0][!grepl('GO_', ONTOLOGY),list(SYMBOL, ID)]
  setkey(nongo)
  nongo <- unique(nongo)
  nongo[,PID:=ID]
  godt <- data.table(Parent=rep(names(func_ann$go_offspring), sapply(func_ann$go_offspring, length)),
                     Child=unname(unlist(func_ann$go_offspring)))
  godt <- rbindlist(list(godt, data.table(Parent=godt[!is.na(Child),unique(Parent)],Child=godt[!is.na(Child),unique(Parent)])))
  godt[is.na(Child),Child:=Parent]
  go <- func_ann$annotations[.(ids),,nomatch = 0][grepl('GO_', ONTOLOGY),list(SYMBOL, ID)]
  setkey(go)
  go <- unique(go)
  toadd <- setdiff(go[,ID],godt[,Child])
  toadddt <- data.table(Parent=toadd,Child=toadd)
  godt <- rbindlist(list(godt, toadddt))
  setkey(go, ID)
  setkey(godt, Child)
  gores <- godt[go,,allow.cartesian = T]
  setnames(gores, c('Parent', 'Child'), c('PID', 'ID'))
  setcolorder(gores, c(identifier, 'ID', 'PID'))
  # setkey(godt, Parent)
  # gores <- rbindlist(pbmclapply(intersect(godt[,unique(Parent)], go[,unique(ID)]), function(i) {
  #   tmp <- go[.(godt[.(i),Child]),,nomatch = 0]
  #   tmp[,PID:=i]
  #   return(unique(tmp))
  # }, mc.cores = mc.cores))
  idlist <- rbindlist(list(nongo, gores))
  return(idlist)
}

jaccard_similarity <- function(m) {
  A <- tcrossprod(m)
  im <- which(A > 0, arr.ind=TRUE, useNames = F)
  b <- rowSums(m)
  Aim <- A[im]
  sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
}
jaccard_distance <- function(m) {
  1 - jaccard_similarity(m)
}

jacdist <- function(identifier, func_ann) {
  tmp <- generateTermAssignment(ids = unique(func_ann$annotations[,eval(as.name(identifier))]), func_ann = func_ann, identifier = identifier)
  tmp <- tmp[!is.na(eval(as.name(identifier)))]
  tmp <- tmp[,list(SYMBOL,PID)]
  tmp <- unique(tmp)
  tmp[,tmp:=1]
  cst <- acast(data = tmp, formula = PID~eval(as.name(identifier)), value.var = 'tmp', fill = 0)
  res <- matrix(c(F,T)[c(cst+1)], nrow = nrow(cst))
  dis <- jaccard_distance(m = res)
  rownames(dis) <- rownames(cst)
  colnames(dis) <- rownames(cst)
  return(dis)
}

jacdist_o <- function(identifier, func_ann) {
  require(pbmcapply)
  tmp <- func_ann$annotations[,c('ID', identifier),with=F][,tmp:=1]
  setkey(tmp)
  tmp <- unique(tmp)
  tmp <- tmp[!is.na(eval(as.name(identifier)))]
  cst <- acast(data = tmp, formula = ID~eval(as.name(identifier)), value.var = 'tmp', fill = 0)
  res <- do.call(cbind, pbmclapply(colnames(cst), function(i) {
    if (i%in%names(func_ann$go_offspring)) {
      as.logical(rowSums(cst[,colnames(cst)%in%func_ann$go_offspring[[i]],drop=F]))
    } else {
      as.logical(cst[,i])
    }
  }, mc.cores = 20))
  dis <- jaccard_distance(m = res)
  rownames(dis) <- rownames(cst)
  colnames(dis) <- rownames(cst)
  return(dis)
}
func_ann <- readRDS('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/preproc/hsa_ann.rds')
dis <- jacdist(identifier = 'SYMBOL', func_ann)
# rownames(dis) <- rownames(cst)
# colnames(dis) <- rownames(cst)
# setkey(func_ann$names, Name)
# dupids <- func_ann$names[.(func_ann$names[,length(unique(ID)),by=Name][V1>1,Name]),ID]
# dis <- dis[!rownames(dis)%in%dupids,!colnames(dis)%in%dupids]
# dis <- as.matrix(dis[!is.na(func_ann$names[.(rownames(dis)),Name]),!is.na(func_ann$names[.(colnames(dis)),Name])])
# setkey(func_ann$names, ID)
# rownames(dis) <- func_ann$names[.(rownames(dis)),Name]
# colnames(dis) <- func_ann$names[.(colnames(dis)),Name]
saveRDS(dis, '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/jaccard_new.rds', compress = T)

parseMsigDB <- function(pathtogmt) {
  gmt <- readLines(pathtogmt)
  gmt <- strsplit(gmt, "\t")
  names <- sapply(gmt, "[", 1)
  gmt <- sapply(gmt, "[", -(1:2))
  names(gmt) <- names
  gmt <- rbindlist(lapply(seq_along(gmt), function(i) data.table(ID=names(gmt)[i],SYMBOL=gmt[[i]])))
  return(gmt)
}

hall <- parseMsigDB("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Dat/MsigDB/20180625/h.all.v6.1.symbols.gmt")
pth <- parseMsigDB("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Dat/MsigDB/20180625/c2.cp.v6.1.symbols.gmt")

msigdb <- rbindlist(list(hall, pth))
func_ann <- list()
func_ann$annotations <- msigdb
func_ann$go_offspring <- list()
dis <- jacdist(identifier = 'SYMBOL', func_ann)
saveRDS(dis, '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/msigdb.rds', compress = T)
