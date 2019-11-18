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

mat1 <- readRDS("../martin/res/20180315_wgcna_allsites_allproteins_incnas_nci60/RDS/mat1.rds")
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
fwrite(termlistdt, '~/phospho/Manuscript/current/supplementTables/new/terms_nci60_noids.txt', sep = '\t', col.names = T, row.names = F)

# pg <- readRDS("Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")

emat <- apply(mat1, 1, function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE) - log10(2)
  x
})
emat <- t(emat)

scores <- gsva(emat, termlist)

smat <- scores$es.obs
# i <- fmatch(rownames(smat), func_ann$names$ID)
# x <- func_ann$names$Name[i]
# rownames(smat) <- x
# 
# grep("cell cycle", x, value = TRUE)
# 
# v <- sort(smat["regulation of cell cycle", ])
# barplot(v, las = 2, col = pg$xref[names(v), "color"])

if (!exists("Res/20180503_gsva/"))
  dir.create("Res/20180503_gsva/")
saveRDS(smat, file = "Res/20180503_gsva/gsva_nci60_trms_new2.RDS")
smat <- readRDS("Res/20180503_gsva/gsva_nci60_trms_new2.RDS")

source("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigures2/R/colorScheme.R")
design <- .cs$nciColor
for (i in seq_along(.cs$tooColor)) {
  design[design==.cs$tooColor[i]] <- names(.cs$tooColor)[i]
}

source("../martin/src/limmafun.R")
mods <- fread("../martin/res/20180315_wgcna_allsites_allproteins_incnas_nci60/CSVs/wgcna_NoImputation.csv")
tl2 <- split(mods$id, mods$moduleColor)
scr <- gsva(emat, tl2) # ,  no.bootstraps=1e4
smt <- scr$es.obs
saveRDS(smt, file = "Res/20180503_gsva/gsva_nci60_mods.RDS")
smt <- readRDS("Res/20180503_gsva/gsva_nci60_mods.RDS")

lim <- limmafun(mat = smat, design = design, p.adjust = "BH", cutoff = 0.05)
edg1 <- parselim(lim)
setnames(edg1, c('Gene set ID', 'Tissue of origin', 'q-value', 'log2-fold-change'))
edg1[`Tissue of origin`=='CN',`Tissue of origin`:='CNS']
fwrite(edg1, '~/phospho/Manuscript/current/supplementTables/new/terms_to_cellline_nci60.txt', sep = '\t', col.names = T, row.names = F)
setnames(edg1, c('source', 'target', 'weight', 'lfc'))
edg1[,target:=.cs$tooColor[target]]
saveRDS(edg1[,list(source,target,weight)], file = "Res/20180503_gsva/edg1_nci60_new.RDS")

lim2 <- limmafun(mat = smt, design = design, p.adjust = "BH", cutoff = 0.05)
edg2 <- parselim(lim2)
setnames(edg2, c('Module', 'Tissue of origin', 'q-value', 'log2-fold-change'))
edg2[`Tissue of origin`=='CN',`Tissue of origin`:='CNS']
fwrite(edg2, '~/phospho/Manuscript/current/supplementTables/new/modules_to_cellline_nci60.txt', sep = '\t', col.names = T, row.names = F)
setnames(edg2, c('source', 'target', 'weight', 'lfc'))
edg2[,target:=.cs$tooColor[target]]
saveRDS(edg2[,list(source,target,weight)], file = "Res/20180503_gsva/edg2_nci60.RDS")
