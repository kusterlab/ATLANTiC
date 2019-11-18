setwd("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen")
library(GSVA)
library(fastmatch)
library(data.table)

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

mat1 <- readRDS("../martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/RDS/mat1.rds")
gname <- sapply(strsplit(rownames(mat1), "_"), "[", 1)
gnl <- split(rownames(mat1), gname)

func_ann <- readRDS('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/preproc/hsa_ann.rds')
terms <- generateTermAssignment(unique(gname), func_ann, 'SYMBOL')
termlist <- split(terms$SYMBOL, terms$PID)
tll <- sapply(termlist, length)
# terms <- func_ann$annotations
# terms <- terms[!is.na(terms$SYMBOL), ]
# terms <- terms[terms$SYMBOL %in% names(gnl), ]
# termlist <- split(terms$SYMBOL, terms$ID)
# tll <- sapply(termlist, length)

# termlist <- termlist[tll > 15 & tll < 500]
termlist <- lapply(termlist, function(x) {
  unique(unlist(gnl[x]))
})

termlistdt <- data.table(`Gene set ID`=rep(names(termlist), sapply(termlist, length)),
                         `Gene name/p-site`=unlist(termlist))
fwrite(termlistdt, '~/phospho/Manuscript/current/supplementTables/new/terms_crc65_noids.txt', sep = '\t', col.names = T, row.names = F)

emat <- apply(mat1, 1, function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE) - log10(2)
  x
})
emat <- t(emat)

scores <- gsva(emat, termlist) # ,  no.bootstraps=1e4

smat <- scores$es.obs
# i <- fmatch(rownames(smat), func_ann$names$ID)
# x <- func_ann$names$Name[i]
# rownames(smat) <- x


if (!exists("Res/20180503_gsva/"))
  dir.create("Res/20180503_gsva/")
saveRDS(smat, file = "Res/20180503_gsva/gsva_crc65_trms_new2.RDS")
smat <- readRDS("Res/20180503_gsva/gsva_crc65_trms_new2.RDS")

traits_crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/traits_crc65.rds')
design <- c('MSI', 'MSS')[traits_crc65[order(traits_crc65$`MSI stable`),"MSI stable",drop=F][colnames(smat),]+1]
names(design) <- colnames(smat)

source("../martin/src/limmafun.R")
mods <- fread("../martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/CSVs/wgcna_NoImputation.csv")
tl2 <- split(mods$id, mods$moduleColor)
scr <- gsva(emat, tl2) # ,  no.bootstraps=1e4
smt <- scr$es.obs
saveRDS(smt, file = "Res/20180503_gsva/gsva_crc65_mods.RDS")
smt <- readRDS("Res/20180503_gsva/gsva_crc65_mods.RDS")

lim <- limmafun(mat = smat, design = design, p.adjust = "BH", cutoff = 0.05)
edg1 <- parselim(lim)
saveRDS(edg1[,list(source,target,weight)], file = "Res/20180503_gsva/edg1_crc65_new.RDS")
setnames(edg1, c('Gene set ID', 'MSI status', 'q-value', 'log2-fold-change'))
edg1[`MSI status`=='MSI',`MSI status`:='+']
edg1[`MSI status`=='MSS',`MSI status`:='-']
fwrite(edg1, '~/phospho/Manuscript/current/supplementTables/new/terms_to_cellline_crc65.txt', sep = '\t', col.names = T, row.names = F)

lim2 <- limmafun(mat = smt, design = design, p.adjust = "BH", cutoff = 0.05)
edg2 <- parselim(lim2)
saveRDS(edg2[,list(source,target,weight)], file = "Res/20180503_gsva/edg2_crc65.RDS")
setnames(edg2, c('Module', 'MSI status', 'q-value', 'log2-fold-change'))
edg2[`MSI status`=='MSI',`MSI status`:='+']
edg2[`MSI status`=='MSS',`MSI status`:='-']
fwrite(edg2, '~/phospho/Manuscript/current/supplementTables/new/modules_to_cellline_crc65.txt', sep = '\t', col.names = T, row.names = F)
