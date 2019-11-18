require(HiveR)
require(data.table)
require(reshape2)
require(scales)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(grid)
require(dendsort)
bgcolor <- alpha('grey90', alpha = 0.3)
ndcolor <- alpha('grey70', alpha = 1)

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

source('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/src/edge2HPD_FIX.R')
source("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigures2/R/colorScheme.R")
alldup <- function (value) { 
  duplicated(value) | duplicated(value, fromLast = TRUE)
}
traits_crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/traits_crc65.rds')
cells <- rownames(traits_crc65[order(traits_crc65$`MSI stable`),"MSI stable",drop=F])

# gsva <- readRDS("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180503_gsva/gsva_crc65.RDS")
# gsva <- gsva[which(!alldup(rownames(gsva))),]
# gsva <- t(scale(t(gsva)))
# # gsva_clust <- dendsort(hclust(as.dist(1-cor(t(gsva), use = 'p')), method = 'ward.D2'), type = 'min')
# # gsva_order <- gsva_clust$labels[gsva_clust$order]
# gsva <- melt(gsva)
# # colnames(gsva) <- c('target', 'source', 'weight')
# colnames(gsva) <- c('source', 'target', 'weight')
# gsva <- gsva[,c('source', 'target', 'weight')]
# gsva <- gsva[gsva$weight>2,]
# gsva <- gsva[!is.na(gsva$source),]
# gsva$weight <- rescale(gsva$weight, to = c(0,1))
# gsva$source <- as.character(gsva$source)
# gsva$target <- as.character(gsva$target)
gsva <- as.data.frame(readRDS("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180503_gsva/edg1_crc65_new.RDS"))

modulemembers <- rbindlist(lapply(list.files('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/CSVs/Modules/NoImputation/', pattern = '.csv', full.names = T), function(i) fread(i, stringsAsFactors = F, integer64 = 'double', header = F)[,Module:=gsub('.csv', '', basename(i))]))
setnames(modulemembers, 'V1', 'id')

terms <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/CSVs/GO/goall_new_rollup_NoImputation.csv')
# setkey(terms, Name)
# terms <- terms[.(terms[,length(unique(ID)),by=Name][V1==1,Name]),]
# terms <- terms[!is.na(Name),]
# terms[grepl('CORUM', ID),ID:=unlist(lapply(strsplit(ID, '_'), function(i) paste(i[1], 'hsa', i[2], sep = '_')))]
# termcount <- terms[,list(`count in term`=length(unique(unlist(strsplit(overlap, ';'))))),by=ID]
termcount <- terms[,list(ID, `count in term`)]
setkey(termcount)
termcount <- unique(termcount)
setkey(termcount, ID)
terms <- data.frame(terms[,list(ID, variable, `q value`)])
colnames(terms) <- c('source', 'target', 'weight')
terms$weight <- rescale(-log10(terms$weight), to = c(0,1))

# termstocellline <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/CSVs/Module_Sample_Correlation_NoImputation.csv')
# termstocellline <- termstocellline[`p value`<=0.05,]
# termstocellline <- data.frame(termstocellline[,list(gsub('^ME', '', Module), Sample, `Correlation`)])
# colnames(termstocellline) <- c('source', 'target', 'weight')
# termstocellline$weight <- rescale(termstocellline$weight, to = c(0,1))
termstocellline <- as.data.frame(readRDS("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180503_gsva/edg2_crc65.RDS"))

modulecount <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/CSVs/wgcna_NoImputation.csv')
modulecount <- modulecount[,length(unique(id)),by=moduleColor]
setkey(modulecount, moduleColor)

moduledist <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/RDS/METreelist.rds')
moduleorder <- moduledist$NoImputation$labels[moduledist$NoImputation$order]
moduleorder <- gsub('ME', '', moduleorder)

tt <- do.call(rbind, list(gsva, terms, termstocellline))
tt <- tt[!is.na(tt$source) & !is.na(tt$target) & !is.na(tt$weight),]
# tt$source <- as.character(tt$source)
# tt$target <- as.character(tt$target)

hpd <- edge2HPD_FIX(tt, axis.cols = alpha('white', alpha = 0), type = '2D', desc = 'Test') # turn edgelist into a workable object
rad <- c(0.3,0.8)
hpd$nodes[hpd$nodes$lab%in%gsva$source, "axis"] = as.integer(2) # Assign all terms to axis 3
hpd$nodes[hpd$nodes$lab%in%gsva$target, "axis"] = as.integer(3) # Assign all cell lines to axis 3
hpd$nodes[hpd$nodes$lab%in%terms$source, "axis"] = as.integer(2) # Assign all terms to axis 2
hpd$nodes[hpd$nodes$lab%in%terms$target, "axis"] = as.integer(1) # Assign all modules to axis 1
hpd$nodes[hpd$nodes$lab%in%termstocellline$source, "axis"] = as.integer(1) # Assign all modules to axis 3
hpd$nodes[hpd$nodes$lab%in%termstocellline$target, "axis"] = as.integer(3) # Assign all cell lines to axis 2
# traits_crc65[order(traits_crc65$`MSI stable`),"MSI stable",drop=F][hpd$nodes[hpd$nodes$axis%in%3, "lab"],]
# hpd$nodes[hpd$nodes$axis%in%3, "color"] <- .cs$colorScheme(scheme = 'c2_2')[traits_crc65[order(traits_crc65$`MSI stable`),"MSI stable",drop=F][hpd$nodes[hpd$nodes$axis%in%3, "lab"],]+1] # Assign node colors for cell lines
hpd$nodes[hpd$nodes$axis%in%3, "color"] <- .cs$colorScheme(scheme = 'c2_1')
# hpd$nodes[hpd$nodes$axis%in%1, "color"] <- hpd$nodes[hpd$nodes$axis%in%1, "lab"] # Assign node colors for modules
hpd$nodes[hpd$nodes$axis%in%1, "color"] <- ndcolor # Assign node colors for modules
hpd$nodes[hpd$nodes$axis%in%2, "color"] <- ndcolor # Assign node color for terms
hpd$nodes[hpd$nodes$axis%in%1, "size"] <- rescale(modulecount[.(hpd$nodes[hpd$nodes$axis%in%1, "lab"]),V1], to = rad) # Assign node size for modules
termc <- termcount[.(hpd$nodes[hpd$nodes$axis%in%2, "lab"]),`count in term`]
termc[is.na(termc)] <- min(termc, na.rm = T)
termc <- rescale(termc, to = rad)
hpd$nodes[hpd$nodes$axis%in%2, "size"] <- termc # Assign node color for terms
hpd$nodes[hpd$nodes$axis%in%3, "size"] <- min(termc, na.rm = T) # Assign node size for cell lines

termdis <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/jaccard_new.rds')
func_ann <- readRDS('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/preproc/hsa_ann.rds')
setkey(func_ann$names, ID)
trms <- hpd$nodes[hpd$nodes$axis%in%2, "lab"]
# IDs <- func_ann$names[.(trms),ID,mult='first']
# termdis <- as.dist(as.matrix(termdis)[rownames(termdis)%in%trms,colnames(termdis)%in%trms])
termdis <- as.dist(as.matrix(termdis)[trms,trms])
term_clust <- dendsort(hclust(termdis, method = 'ward.D2'), type = 'min')
term_order <- term_clust$labels[term_clust$order]

rad <- c(0.01,0.3)
# hpd$nodes[hpd$nodes$axis%in%3, "radius"] <- rescale(match(hpd$nodes[hpd$nodes$axis%in%3, "lab"], cells), to = rad) # Assign node colors for cell lines
hpd$nodes[hpd$nodes$axis%in%3, "radius"] <- rescale(c(0:2), to = rad)[2:3] # Assign node colors for cell lines
hpd$nodes[hpd$nodes$axis%in%1, "radius"] <- rescale(match(hpd$nodes[hpd$nodes$axis%in%1, "lab"], moduleorder), to = rad) # Assign node colors for modules
hpd$nodes[hpd$nodes$axis%in%2, "radius"] <- rescale(match(hpd$nodes[hpd$nodes$axis%in%2, "lab"], term_order), to = rad) # Assign node color for terms
# ln <- length(hpd$nodes[hpd$nodes$axis%in%2, "radius"][is.na(hpd$nodes[hpd$nodes$axis%in%2, "radius"])])
# hpd$nodes[hpd$nodes$axis%in%2, "radius"][is.na(hpd$nodes[hpd$nodes$axis%in%2, "radius"])] <- seq(max(hpd$nodes[hpd$nodes$axis%in%2, "radius"], na.rm = T)+1, max(hpd$nodes[hpd$nodes$axis%in%2, "radius"], na.rm = T)+ln, by = 1)
# hpd$nodes[hpd$nodes$axis%in%2, "radius"] <- rescale(hpd$nodes[hpd$nodes$axis%in%2, "radius"], to = rad)

circlefinder <- function(hpd) {
  require(colorspace)
  nodes <- hpd$nodes[,c('id','axis','color')]
  edges <- hpd$edges[,c('id1','id2','weight')]
  setDT(nodes)
  setDT(edges)
  setkey(nodes, id)
  cc <- edges[,axis1:=nodes[.(id1),axis]]
  cc[,axis2:=nodes[.(id2),axis]]
  tbl1 <- cc[axis1==1 & axis2 == 3,list(id1, id2)]
  setkey(tbl1, id2)
  tbl2 <- cc[axis1==2 & axis2 == 3,list(id1, id2)]
  setkey(tbl2, id2)
  tbl <- merge(tbl1, tbl2, by = 'id2', allow.cartesian = T)
  # setkey(tbl, id1.y, id1.x)
  # setkey(cc, id1, id2)
  test <- merge(tbl,cc, by.x = c("id1.y", "id1.x"), by.y = c("id1", "id2"), all = F)
  setkey(test, id2)
  setkey(nodes, id)
  test <- test[nodes,,nomatch = 0]
  res1 <- test[,list(id1=id1.y, id2=id1.x, color)]
  res2 <- test[,list(id1=id1.y, id2=id2, color)]
  res3 <- test[,list(id1=id1.x, id2=id2, color)]
  res <- rbindlist(list(res1,res2,res3))
  setkey(res)
  res <- unique(res)
  setnames(res, 'color', 'circlecolor')
  # dupids <- res[,length(unique(color)),by=list(id1,id2)][V1>1,list(id1,id2)]
  # res[paste0(id1, '_', id2)%in%dupids[,paste0(id1, '_', id2)],color:='pink']
  # setkey(res)
  # res <- unique(res)
  return(res)
  # tst <- data.table(id1 = c(10000,10001,10002),
  #                   id2 = c(10001,10002,10000),
  #                   axis1 = c(2, 3, 1),
  #                   axis2 = c(1, 2, 3))
  # cc <- rbindlist(list(cc, tst))
  # setkey(cc, id1)
  # ids <- cc[axis1==2 & axis2==1,unique(id1)]
  # for (i in ids) {
  #   i2 <- cc[.(i),id2,nomatch=0]
  #   i3 <- cc[.(unique(i2)),id2,nomatch=0]
  #   i4 <- cc[.(unique(i3)),id2,nomatch=0]
  #   if (length(intersect(i, unique(i4)))==1) {
  #     cc[id2%in%i & id1%in%i3,cand:=1]
  #     tmp <- cc[id2%in%i & id1%in%i3,unique(id1)]
  #     cc[id2%in%tmp & id1%in%i2,cand:=1]
  #     tmp <- cc[id2%in%tmp & id1%in%i2,unique(id1)]
  #     cc[id2%in%tmp & id1%in%i,cand:=1]
  #   }
  # }
  # return(cc[cand==1,])
}
cc <- circlefinder(hpd)

# color_lookup = c("Protein" = "Greens", "Regulator" = "Reds", "Substrate" = "Blues") # assign each type a unique color palette (will be used with ColorBrewer below)
# color_counts = c("Protein" = 4, "Regulator" = 3, "Substrate" = 3) # these numbers are used in ColorBrewer to determine the number of colors in each category
# 
# # Helper function to color the edges according to the origin
# determine_line_color <- function(id1)
# {
#   source = hpd$nodes[id1, "lab"]
#   
#   type = str_sub(source, 1, -3) # extract the "Protein", "Regulator" or "Substrate" part
#   number = as.numeric(str_sub(source, -1, -1)) # extraxt the "id"
#   
#   color = brewer.pal(color_counts[type], color_lookup[type])[number]
#   
#   return(color)
# }
# 
# hpd$edges$color <- sapply(hpd$edges$id1, determine_line_color) # assign colors to the lines based on the source
hpd$edges$color <- bgcolor
hpd$edges <- merge(hpd$edges, cc, by.x = c('id1', 'id2'), by.y = c('id1', 'id2'), all = T)
hpd$edges[!is.na(hpd$edges$circlecolor),'color'] <- hpd$edges[!is.na(hpd$edges$circlecolor),'circlecolor']
hpd$edges$circlecolor <- NULL
# tmp <- paste0(hpd$edges$id1, '_', hpd$edges$id2)
# tmp2 <- paste0(cc$id1, '_', cc$id2)
# cc[,ids:=tmp2]
# setkey(cc, ids)
# ids <- which(tmp%in%tmp2)
# idx <- tmp[ids]
# # hpd$edges[ids,'color'] <- alpha('green', alpha = 0.3)
# hpd$edges[ids,'color'] <- cc[.(idx),color]
hpd$edges <- hpd$edges[c(which(hpd$edges$color==bgcolor),which(hpd$edges$color!=bgcolor)),]
hpd$edges$weight <- 0.1
# identical(bckp,hpd$edges)
# hpd$edges[hpd$edges$id1==1 & hpd$edges$id2==2720,'color'] <- 'green'
# hpd$edges[hpd$edges$id1==2720 & hpd$edges$id2%in%c(2849,2931),'color'] <- 'green'
# hpd$edges[hpd$edges$id1%in%c(2849,2931) & hpd$edges$id2%in%c(1),'color'] <- 'green'

# redhpd <- hpd
# set.seed(1)
# ids <- which(hpd$edges$color=='green')
# redhpd$edges <- hpd$edges[c(sample(1:nrow(hpd$edges[-ids,]), size = 1000, replace = F), ids),]
# redhpd$edges <- hpd$edges[c(sample(1:nrow(hpd$edges), size = 1000, replace = F)),]

# Create the hive plot
dir.create(path = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180506_hive_crc65', showWarnings = F)
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180506_hive_crc65/test.png', res = 300)
# pdf(file = paste0('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180506_hive_crc65/test.pdf'), width = 5.7/2.539998, height = 5.7/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular', useDingbats = F, encoding = 'Greek')
# plotHive(tmp, method = 'abs',
#          axLabs = c("Modules", "Terms", "MSI"), 
#          bkgnd = "white", anNodes = 'nodes.txt',
#          axLab.gpar = gpar(fontsize = 10, col = 'black', lwd = 10),ch = 0, axLab.pos = 0.06, rot = c(0, -60, 60, c(rep(0, nrow(fread('nodes.txt'))))))
# # plotHive(hpd,
# #          axLabs = c("Modules", "Terms", "Cell lines"), 
# #          bkgnd = "white", 
# #          axLab.gpar = gpar(col = "#bbbbbb"),ch = 0)
# dev.off()

genLabels <- function(hpd, labs, maxradius, offset, func_ann, filename = 'nodes.txt', dclust = T) {
  require(mclust)
  require(data.table)
  res <- hpd$nodes[hpd$nodes$lab%in%labs,]
  setDT(res)
  setkey(func_ann$names, ID)
  setkey(res, axis, radius)
  if (dclust) {
    try(res[axis%in%res[,.N,by=axis][N>1,axis],grp:=paste0(color, '_', densityMclust(radius, modelNames = 'E')$classification),by=list(axis,color)], silent = T)
    res[axis%in%res[,.N,by=axis][N==1,axis],grp:=1]
  } else {
    res[,grp:=1:.N,by=list(axis)]
  }
  res[,name:=lab]
  res[axis==2,name:=func_ann$names[.(lab),Name]]
  ltbl <- res[,list(node.lab=lab[ceiling(length(unique(radius))/2)],
                    node.text=name[ceiling(length(unique(radius))/2)]),by=list(axis,grp)]
  setkey(res, lab)
  ltbl[,node.text:=1:.N]
  ltbl[,grp:=NULL]
  ltbl[,angle:=1]
  ltbl[,radius:=(maxradius-res[.(node.lab),radius])+0.05]
  ltbl[,offset:=offset]
  ltbl[,hjust:=0.5]
  ltbl[,vjust:=0.5]
  ltbl[axis==1,angle:=360]
  ltbl[axis==2,angle:=120]
  ltbl[axis==3,angle:=240]
  ltbl[,axis:=NULL]
  # require(tm)
  # require(topicmodels)
  # docs <- Corpus(VectorSource(res[axis==2,lab]))
  # #remove potentially problematic symbols
  # toSpace <- content_transformer(function(x, pattern) { return (gsub(pattern, " ", x))})
  # docs <- tm_map(docs, toSpace, "-")
  # docs <- tm_map(docs, toSpace, "’")
  # docs <- tm_map(docs, toSpace, "‘")
  # 
  # docs <- tm_map(docs,content_transformer(tolower))
  # docs <- tm_map(docs, removePunctuation)
  # dtm <- DocumentTermMatrix(docs)
  # bla <- LDA(x = dtm, k = max(res[axis==2,grp]))
  # as.matrix(terms(bla,10))
  # as.matrix(topics(bla))
  
  fwrite(ltbl, file = filename)
}

# require(pbmcapply)
cols <- c(as.list(unique(hpd$edges$color)[!unique(hpd$edges$color)%in%bgcolor]),list(unique(hpd$edges$color)[!unique(hpd$edges$color)%in%bgcolor]))
anns <- lapply(rev(cols), function(i) {
  tmp <- hpd
  tmp$edges[!tmp$edges$color%in%i,'color'] <- bgcolor
  tmp$edges <- tmp$edges[c(which(tmp$edges$color%in%bgcolor),which(!tmp$edges$color%in%bgcolor)),]
  tmp$nodes[!tmp$nodes$color%in%i,'color'] <- ndcolor
  for (j in i) {
    tmp$nodes[tmp$nodes$id%in%unique(unlist(tmp$edges[tmp$edges$color==j,c('id1','id2')])),'color'] <- j
  }
  tmp$nodes <- tmp$nodes[c(which(tmp$nodes$color%in%ndcolor),which(!tmp$nodes$color%in%ndcolor)),]
  # tmp$edges <- tmp$edges[tmp$edges$color!=bgcolor,]
  
  labs <- tmp$nodes$lab[!tmp$nodes$color%in%ndcolor] #  & tmp$nodes$axis==2
  
  # filter by Name
  bla2 <- merge(func_ann$names[.(labs)], hpd$nodes, by.x = 'ID', by.y = 'lab', all.x = T)
  setkey(bla2, color, radius)
  bla3 <- bla2[grepl('mis|mmr|msh|mll|muts', ignore.case = T, Name) | grepl('rna|transc', ignore.case = T, Name),ID]
  bla4 <- bla2[is.na(Name),ID]
  bla2 <- merge(func_ann$names[.(c(bla3,bla4))], hpd$nodes, by.x = 'ID', by.y = 'lab', all.x = T)
  setkey(bla2, id)
  bla1 <- bla2[.(unique(tmp$edges[tmp$edges$id1%in%bla2[axis==2,id] & tmp$edges$id2%in%bla2[axis%in%c(1,3),id],'id2'])),ID]
  bla4 <- intersect(bla1, bla4)
  labs <- c(bla3,bla4)[match(labs, c(bla3,bla4), nomatch = 0)]
  
  nds <- file.path('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180506_hive_crc65', paste0(paste0(i, collapse = '_'), '.txt'))
  genLabels(hpd = tmp, labs = labs, maxradius = 0.3, offset = 0.03, filename = nds, func_ann = func_ann, dclust = F)
  
  # setEPS()
  pdf(file = paste0('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180506_hive_crc65/test_', paste0(i, collapse = '_'), '_new3.pdf'), width = 5.7/2.539998, height = 5.7/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular', useDingbats = F, encoding = 'Greek')
  # gpar <- gpar(mar = c(0,0,0,0),
  #      oma = c(0,0,0,0),
  #      lwd = 0.01)
  plotHive(tmp, method = 'abs',
           axLabs = c("2) Groups of abundance-\nrelated Proteins/p-sites", "1) Groups of functionally\nrelated Proteins/p-sites", "3) Groups of cell lines\nwith the same MSI status"), 
           bkgnd = "white", anNodes = NULL,
           axLab.gpar = gpar(fontsize = 10, col = 'black', lwd = 10),ch = 0, axLab.pos = 0.06, rot = c(0, -60, 60, c(rep(0, nrow(fread(nds))))))
  dev.off()
  bla <- fread(nds)
  blax <- merge(func_ann$names[.(bla$node.lab),,nomatch=0], bla[,list(node.lab, radius)], by.x = 'ID', by.y = 'node.lab')[order(radius),]
  blay <- merge(bla, hpd$nodes, by.x = 'node.lab', by.y = 'lab')
  setkey(blay, axis, radius.y)
  blay[,term:=func_ann$names[.(node.lab),Name]]
  blay[,list(term,node.lab,node.text,radius.y,radius.x)]
  return(blay[,list(term,node.lab,node.text,radius.y,radius.x)])
})

blub <- rbindlist(lapply(cols[1:2], function(i) {
  cat(i, '\n')
  tmp <- hpd
  # tmp$edges[!tmp$edges$color%in%i,'color'] <- bgcolor # color all edges not of color i in the background color
  # tmp$edges <- tmp$edges[c(which(!tmp$edges$color%in%bgcolor)),]
  tmp$edges <- tmp$edges[tmp$edges$color%in%i,] # select edges of color i
  if (nrow(tmp$edges)==0)
    return(data.table(lab1=character(0),Name=character(0),lab2=character(0)))
  # tmp$nodes[!tmp$nodes$color%in%i,'color'] <- ndcolor
  tmp$nodes <- tmp$nodes[tmp$nodes$id%in%unique(unlist(tmp$edges[,c('id1','id2')])),]
  tmp$nodes$color <- i
  # for (j in i) {
  #   tmp$nodes[tmp$nodes$id%in%unique(unlist(tmp$edges[tmp$edges$color==j,c('id1','id2')])),] <- j
  # }
  # tmp$nodes <- tmp$nodes[c(which(tmp$nodes$color%in%ndcolor),which(!tmp$nodes$color%in%ndcolor)),]
  # tmp$edges <- tmp$edges[tmp$edges$color!=bgcolor,]
  
  setDT(tmp$nodes)
  setDT(tmp$edges)
  
  setkey(tmp$nodes, id)
  setkey(tmp$edges, id1, id2)
  
  # add axis
  res <- tmp$edges[tmp$nodes,on=.(id1==id),nomatch = 0][,list(id1, axis1=axis, color1=color,lab1=lab, id2, axis2=NA, color2=NA, lab2=NA)]
  res <- res[tmp$nodes,on=.(id2==id),nomatch = 0][,list(id1, axis1, color1, lab1, id2, axis2=axis, color2=color,lab2=lab)]
  
  res <- res[(axis1==1 & axis2 == 2) | (axis1==2 & axis2 == 1),]
  setkey(res, lab1)
  res <- res[func_ann$names,,nomatch=0][,list(lab1,Name,lab2)]
  setkey(res, lab2)
  # tisid <- tmp$nodes[lab==i,id]
  # tmp$edges[,list(id1, id2, axis, color)][id1==tisid | id2==tisid,unlist(list(id1,id2))[!unlist(list(id1,id2))%in%tisid]]
  # 
  # 
  # res <- tmp$edges[,list(id1, id2, color)][tmp$nodes[,list(id, lab)],on=.(id==id),nomatch = 0]
  # setkey(res, id1, id2)
  # 
  # res <- res[res,on=.(id1==id2),nomatch=0][,list(id1,id2,color,cluster=lab,term=i.lab)]
  
  # labs <- tmp$nodes$lab
  # bla <- data.table(tmp$nodes)
  # res <- func_ann$names[ID%in%bla[color==i & axis==2,lab],]
  labs <- c('+', '-')
  names(labs) <- c('MSI', 'MSS')
  res[,MSI:=labs[hpd$nodes[hpd$nodes$axis%in%3,"lab"][hpd$nodes[hpd$nodes$axis%in%3,'color']==i]]]
  return(res)
}))

setkey(blub, lab2)
setkey(modulemembers, Module)

res <- blub[modulemembers,,allow.cartesian=T,nomatch=0]
setnames(res, c('Gene set ID', 'Gene set name', 'Module', 'MSI status', 'Gene name/p-site'))
fwrite(res, '~/phospho/Chen/Res/20180523_moreSuppTables/WGCNA_new_annotations_CRC65_new.txt', col.names = T, row.names = F, sep = '\t')

terms <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/CSVs/GO/goall_new_rollup_NoImputation.csv')
terms <- terms[,list(ID, variable, `q value`)]
setnames(terms, c('source', 'target', 'weight'))

# lapply(anns, setkey, node.lab)
# terms[,name:=anns[[3]][.(source),term]]
# mods <- terms[name%in%anns[[3]]$term & target%in%anns[[3]]$node.lab,][!is.na(name),unique(target)] #  & (grepl('mis|mmr|msh|mll|muts', ignore.case = T, name) | grepl('rna|transc', ignore.case = T, name))
# trms <- terms[name%in%anns[[3]]$term & target%in%anns[[3]]$node.lab,][!is.na(name),unique(source)] #  & (grepl('mis|mmr|msh|mll|muts', ignore.case = T, name) | grepl('rna|transc', ignore.case = T, name))
# potgenes <- func_ann$annotations[ID%in%trms,][order(ID),unique(SYMBOL)]

setkey(func_ann$names, ID)
terms[,name:=func_ann$names[.(source),Name]]
mods <- hpd$nodes[hpd$nodes$id%in%unique(unlist(hpd$edges[hpd$edges$color%in%cols[[3]],c('id1','id2')])) & hpd$nodes$axis%in%1,'lab']
tr <- hpd$nodes[hpd$nodes$id%in%unique(unlist(hpd$edges[hpd$edges$color%in%cols[[3]],c('id1','id2')])) & hpd$nodes$axis%in%2,'lab']
trms <- terms[source%in%tr & target%in%mods,][!is.na(name),unique(source)] #  & (grepl('mis|mmr|msh|mll|muts', ignore.case = T, name) | grepl('rna|transc', ignore.case = T, name))
# Alternative:
trms <- res[,unique(`Gene set ID`)]
# All genes which ARE part of terms that are enriched in coloured correlation clusters
terms_crc65_noids <- fread('~/phospho/Manuscript/current/supplementTables/new/terms_crc65_noids.txt', stringsAsFactors = F, integer64 = 'double')
potgenes <- terms_crc65_noids[`Gene set ID`%in%trms,unique(sapply(strsplit(`Gene name/p-site`, '_'), '[', 1))]
# potgenes <- func_ann$annotations[ID%in%trms,][order(ID),unique(SYMBOL)]

modulecount <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/CSVs/wgcna_NoImputation.csv')
setkey(modulecount, moduleColor)
actgenes <- modulecount[mods,][as.logical(rowSums(sapply(potgenes, function(i) grepl(i, id)), na.rm = T)),]
setkey(actgenes, id)
# The genes present in coloured correlation clusters
gns <- modulecount[mods,id]
actgenes2 <- modulecount[mods,][id%in%potgenes,]
setkey(actgenes2, id)

# scr <- readRDS("../Chen/Res/20180503_gsva/gsva_crc65_mods.RDS")
msi <- traits_crc65[,"MSI stable"]
names(msi) <- rownames(traits_crc65)
design <- c('MSI', 'MSS')[msi+1]
names(design) <- names(msi)
mat1 <- readRDS("~/phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/RDS/mat1.rds")
design <- design[names(design)%in%colnames(mat1)]

source('~/phospho/martin/src/limmafun.R')
lim <- limmafun(mat = log2(10^(mat1[gns,])), design = design, p.adjust = "BH", cutoff = 1)
lim$tab[,plot(MSS, -log10(adj.P.Val), xlim = c(-3,0), ylim = c(0,6), col = 'white')]
lim$tab[,text(MSS, -log10(adj.P.Val), labels = rn)]

col <- .cs$colorScheme("c2_1")
source("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigures2/R/colorScheme.R")
setEPS()
postscript("~/phospho/Chen/pubFigure3/figures/newfigure2_A_2_new2.eps", width = 5.8 /2.539998, height = 5.8/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular')
par(mar = c(2.5, 3.1, 2.2, 1.2))

plot(lim$tab$MSS, y = -log10(lim$tab$adj.P.Val), col = "white", cex = 1, axes = FALSE, pch = 19, xlab = "", ylab = "")
box()
# abline(h = c(20, 40, 60, 80, 100), col = "gray", lty = 3)
points(lim$tab$MSS[lim$tab$adj.P.Val>0.05], y = -log10(lim$tab$adj.P.Val)[lim$tab$adj.P.Val>0.05], col = 'grey90', cex = 0.5, pch = 19, lwd = 0)
lim$tab[!rn%in%potgenes & MSS<0 & adj.P.Val<0.05,points(MSS, -log10(adj.P.Val), col = col['TUMblue'], cex = 0.5, pch = 19, lwd = 0)]
lim$tab[!rn%in%potgenes & MSS>0 & adj.P.Val<0.05,points(MSS, -log10(adj.P.Val), col = col['TUMred'], cex = 0.5, pch = 19, lwd = 0)]
lim$tab[rn%in%potgenes & MSS<0 & adj.P.Val<0.05,points(MSS, -log10(adj.P.Val), col = al2hex(alpha(col['TUMblue'], 1)), cex = 0.5, pch = 19, lwd = 0)]
lim$tab[rn%in%potgenes & MSS>0 & adj.P.Val<0.05,points(MSS, -log10(adj.P.Val), col = al2hex(alpha(col['TUMred'], 1)), cex = 0.5, pch = 19, lwd = 0)]

# points(dat$effectsize, y = dat$freq, col = "white", cex = 0.85, pch = "_")

# axis(1, at = seq(-5,5,2), line = 0, cex.axis = 1, padj = 0.5, tck=-0.02, mgp = c(3,0,0))
axis(1, at = seq(-6,6,2), line = 0, cex.axis = 1, padj = 0.5, tck=-0.02, mgp = c(3,0,0))
axis(2, line = 0, cex.axis = 1, padj = 0.5, tck=-0.02, mgp = c(3,0.5,0), las = 2)

# mtext(side = 2, text = c(20, 40, 60, 80, 100), line = 0.2, at = c(20, 40, 60, 80, 100), las = 2, cex = 0.7)
mtext(side = 1, text = "log2-fold-change", cex = 1, line = 1.3)
mtext(side = 2, text = "-log10(q-value)", cex = 1, line = 2)
mtext(side = 3, text = "Proteins/p-sites differentially\nabundant in MSI+/- cells", cex = 1, line = 0.2)
legend('topleft', bty = 'n', title = 'High in MSI:', legend = c('+', '-'), col = col[c('TUMblue','TUMred')], pch = 19, xjust = 0, title.adj = 0.2, ncol = 2)
lines(c(0.8,0.8), c(3.3,17))
lines(c(0.8,5), c(3.3,3.3))

dev.off()

plot(lim$tab$MSS[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>2], y = -log10(lim$tab$adj.P.Val)[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>2], col = "white", cex = 1, axes = FALSE, pch = 19, xlab = "", ylab = "", xlim = c(0.75,4.05))
box()

lim$tab[!rn%in%potgenes & MSS>0.8 & -log10(adj.P.Val)>2,points(MSS, -log10(adj.P.Val), col = col['TUMred'], cex = 0.5, pch = 19, lwd = 0)]
lim$tab[rn%in%potgenes & MSS>0.8 & -log10(adj.P.Val)>2,points(MSS, -log10(adj.P.Val), col = al2hex(alpha(col['TUMred'], 0.5)), cex = 0.5, pch = 19, lwd = 0)]
# points(lim$tab$MSS[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>3.5], y = -log10(lim$tab$adj.P.Val)[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>3.5], col = col['TUMred'], cex = 0.5, pch = 19, lwd = 0)
text(lim$tab$MSS[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>2], y = -log10(lim$tab$adj.P.Val)[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>2]+0.1, labels = lim$tab$rn[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>2], cex = 0.5, pch = 19)


setEPS()
postscript("~/phospho/Chen/pubFigure3/figures/newfigure2_A_3_new2_annot.eps", width = 5.8 /2.539998, height = 5.8/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular')
par(mar = c(2.5, 3.1, 2.2, 1.2))

plot(lim$tab$MSS[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>3.3], y = -log10(lim$tab$adj.P.Val)[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>3.3], col = "white", cex = 1, axes = FALSE, pch = 19, xlab = "", ylab = "", xlim = c(0.75,4.05))
box()

lim$tab[!rn%in%potgenes & MSS>0.8 & -log10(adj.P.Val)>3.3,points(MSS, -log10(adj.P.Val), col = col['TUMred'], cex = 0.5, pch = 19, lwd = 0)]
lim$tab[rn%in%potgenes & MSS>0.8 & -log10(adj.P.Val)>3.3,points(MSS, -log10(adj.P.Val), col = al2hex(alpha(col['TUMred'], 0.5)), cex = 0.5, pch = 19, lwd = 0)]
# points(lim$tab$MSS[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>3.5], y = -log10(lim$tab$adj.P.Val)[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>3.5], col = col['TUMred'], cex = 0.5, pch = 19, lwd = 0)
text(lim$tab$MSS[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>3.3], y = -log10(lim$tab$adj.P.Val)[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>3.3]+0.1, labels = lim$tab$rn[lim$tab$MSS>0.8 & -log10(lim$tab$adj.P.Val)>3.3], cex = 0.5, pch = 19)

# axis(1, at = seq(-1,5,2), line = 0, cex.axis = 1, padj = 0.5, tck=-0.02, mgp = c(3,0,0))
# axis(1, at = seq(-2,6,2), line = 0, cex.axis = 1, padj = 0.5, tck=-0.02, mgp = c(3,0,0))
# axis(2, line = 0, cex.axis = 1, padj = 0.5, tck=-0.02, mgp = c(3,0.5,0), las = 2)

# mtext(side = 2, text = c(20, 40, 60, 80, 100), line = 0.2, at = c(20, 40, 60, 80, 100), las = 2, cex = 0.7)
# mtext(side = 1, text = "log2-fold-change", cex = 1, line = 1.3)
# mtext(side = 2, text = "-log10(q-value)", cex = 1, line = 2)
# mtext(side = 3, text = "Proteins/p-sites differentially\nabundant in MSI +/- cells", cex = 1, line = 0.2)
# legend('topleft', bty = 'n', title = 'High in MSI', legend = c('+', '-'), col = col[c('TUMblue','TUMred')], pch = 19, xjust = 0, title.adj = 0.2, ncol = 2)
dev.off()

blub <- terms[source%in%trms,list(source,name,target)]
setkey(blub)
blub <- unique(blub)
gsva_dt <- setDT(gsva)
setkey(gsva_dt, source)
setkey(blub, source)
blub[gsva_dt,,nomatch=0]
blub <- blub[gsva_dt,,nomatch=0]
# blub[,length(unique(i.target)),by=source]
bla <- blub[,list(name, i.target)]
setkey(bla)
bla <- unique(bla)
View(bla[order(i.target)])


mmrs <- func_ann$names[.(hpd$nodes$lab[hpd$nodes$id%in%unique(cc$id1,cc$id2)]),,nomatch=0][grepl('mis|mmr|msh|mll|muts', ignore.case = T, Name),]
fid <- hpd$edges$id1[which(hpd$edges$id1%in%hpd$nodes$id[which(hpd$nodes$lab%in%mmrs[,ID])])]
mmedge <- hpd$edges[hpd$edges$id1%in%fid,]
hpd$nodes[hpd$nodes$id%in%mmedge$id2,]
merge(mmrs, hpd$nodes[hpd$nodes$id%in%fid,], by.x = 'ID', by.y = 'lab', all = T)[order(radius),]
