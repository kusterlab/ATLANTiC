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

source('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/src/edge2HPD_FIX.R')
source("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigures2/R/colorScheme.R")
alldup <- function (value) { 
  duplicated(value) | duplicated(value, fromLast = TRUE)
}
# cells <- names(sort(.cs$nciColor))
cells <- unique(sort(.cs$nciColor))

# gsva <- readRDS("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180503_gsva/gsva_nci60.RDS")
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
gsva <- as.data.frame(readRDS("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180503_gsva/edg1_nci60_new.RDS"))
require(scales)
gsva$target <- col2hcl(gsva$target)

terms <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_nci60/CSVs/GO/goall_new_rollup_NoImputation.csv')
# setkey(terms, Name)
# terms <- terms[.(terms[,length(unique(ID)),by=Name][V1==1,Name]),]
# terms <- terms[!is.na(Name),]
# terms[grepl('CORUM', ID),ID:=unlist(lapply(strsplit(ID, '_'), function(i) paste(i[1], 'hsa', i[2], sep = '_')))]
termcount <- terms[,list(ID, `count in term`)]
setkey(termcount)
termcount <- unique(termcount)
setkey(termcount, ID)
terms <- data.frame(terms[,list(ID, variable, `q value`)])
colnames(terms) <- c('source', 'target', 'weight')
terms$weight <- rescale(-log10(terms$weight), to = c(0,1))

# termstocellline <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_nci60/CSVs/Module_Sample_Correlation_NoImputation.csv')
# termstocellline <- termstocellline[`p value`<=0.05,]
# termstocellline <- data.frame(termstocellline[,list(gsub('^ME', '', Module), Sample, `Correlation`)])
# colnames(termstocellline) <- c('source', 'target', 'weight')
# termstocellline$weight <- rescale(termstocellline$weight, to = c(0,1))
termstocellline <- as.data.frame(readRDS("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180503_gsva/edg2_nci60.RDS"))
termstocellline$target <- col2hcl(termstocellline$target)

modulemembers <- rbindlist(lapply(list.files('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_nci60/CSVs/Modules/NoImputation/', pattern = '.csv', full.names = T), function(i) fread(i, stringsAsFactors = F, integer64 = 'double', header = F)[,Module:=gsub('.csv', '', basename(i))]))
setnames(modulemembers, 'V1', 'id')

modulecount <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_nci60/CSVs/wgcna_NoImputation.csv')
modulecount <- modulecount[,length(unique(id)),by=moduleColor]
setkey(modulecount, moduleColor)

moduledist <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_nci60/RDS/METreelist.rds')
moduleorder <- moduledist$NoImputation$labels[moduledist$NoImputation$order]
moduleorder <- gsub('ME', '', moduleorder)

tt <- do.call(rbind, list(gsva, terms, termstocellline))
tt <- tt[!is.na(tt$source) & !is.na(tt$target) & !is.na(tt$weight),]
# tt$source <- as.character(tt$source)
# tt$target <- as.character(tt$target)

hpd <- edge2HPD_FIX(tt, axis.cols = alpha('white', alpha = 0), type = '2D', desc = 'Test') # turn edgelist into a workable object
rad <- c(0.3,0.8)
addnode <- data.frame(as.integer(max(hpd$nodes$id)+seq(1,length(col2hcl(cells)[!col2hcl(cells)%in%hpd$nodes$lab]), 1)), col2hcl(cells)[!col2hcl(cells)%in%hpd$nodes$lab], 1L, 1, 1, 'black', stringsAsFactors = F)
dimnames(addnode) <- list(row.names = 1:nrow(addnode), col.names = colnames(hpd$nodes))
hpd$nodes <- rbind(hpd$nodes, addnode)
hpd$nodes[hpd$nodes$lab%in%addnode$lab, "axis"] = as.integer(3) # Assign all cell lines to axis 2
# addedge <- data.frame(rep(hpd$edges[hpd$edges$id1%in%hpd$nodes[hpd$nodes$axis!=3,'id'],'id1'][1],nrow(addnode)), addnode$id, 0.1, 'white')
# dimnames(addedge) <- list(row.names = 1:nrow(addedge), col.names = colnames(hpd$edges))
# hpd$edges <- rbind(hpd$edges, addedge)
chkHPD(hpd)
hpd$nodes[hpd$nodes$lab%in%gsva$source, "axis"] = as.integer(2) # Assign all terms to axis 3
hpd$nodes[hpd$nodes$lab%in%gsva$target, "axis"] = as.integer(3) # Assign all cell lines to axis 3
hpd$nodes[hpd$nodes$lab%in%terms$source, "axis"] = as.integer(2) # Assign all terms to axis 2
hpd$nodes[hpd$nodes$lab%in%terms$target, "axis"] = as.integer(1) # Assign all modules to axis 1
hpd$nodes[hpd$nodes$lab%in%termstocellline$source, "axis"] = as.integer(1) # Assign all modules to axis 3
hpd$nodes[hpd$nodes$lab%in%termstocellline$target, "axis"] = as.integer(3) # Assign all cell lines to axis 2
# hpd$nodes[hpd$nodes$axis%in%3, "color"] <- .cs$nciColor[hpd$nodes[hpd$nodes$axis%in%3, "lab"]] # Assign node colors for cell lines
hpd$nodes[hpd$nodes$axis%in%3, "color"] <- hpd$nodes[hpd$nodes$axis%in%3, "lab"] # Assign node colors for cell lines
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
termdis <- as.dist(as.matrix(termdis)[trms,trms])
term_clust <- dendsort(hclust(termdis, method = 'ward.D2'), type = 'min')
term_order <- term_clust$labels[term_clust$order]

rad <- c(0.01,0.3)
hpd$nodes[hpd$nodes$axis%in%3, "radius"] <- rescale(match(hpd$nodes[hpd$nodes$axis%in%3, "lab"], col2hcl(cells)), to = rad) # Assign node colors for cell lines
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
dir.create(path = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180506_hive', showWarnings = F)
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180506_hive/test.png', res = 300)
# plotHive(hpd,
#          axLabs = c("Modules", "Terms", "Cell lines"), 
#          bkgnd = "white", 
#          axLab.gpar = gpar(col = "#bbbbbb"),ch = 0)
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
cols <- c(as.list(unique(hpd$nodes$color)[!unique(hpd$nodes$color)%in%ndcolor]),list(unique(hpd$nodes$color)[!unique(hpd$nodes$color)%in%ndcolor]))
anns <- lapply(cols, function(i) {
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
  nds <- file.path('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180506_hive', paste0(paste0(i, collapse = '_'), '.txt'))
  genLabels(hpd = tmp, labs = labs, maxradius = 0.3, offset = 0.03, func_ann = func_ann, filename = nds, dclust = T)
  
  # setEPS()
  pdf(file = paste0('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180506_hive/test2_', paste0(i, collapse = '_'), 'new3.pdf'), width = 5.7/2.539998, height = 5.7/2.539998, pointsize = 8, family = 'TUM Neue Helvetica 55 Regular', useDingbats = F, encoding = 'Greek')
  # gpar <- gpar(mar = c(0,0,0,0),
  #      oma = c(0,0,0,0),
  #      lwd = 0.01)
  plotHive(tmp, method = 'abs',
           axLabs = c("2) Groups of abundance-\nrelated Proteins/p-sites", "1) Groups of functionally\nrelated Proteins/p-sites", "3) Groups of cell lines from\nthe same tissue of origin"), 
           # axLabs = c("Abundance-related\nProteins/p-sites", "Functionally related\nProteins/p-sites", "Cell lines from the\nsame tissue of origin"), 
           bkgnd = "white", anNodes = NULL,
           axLab.gpar = gpar(fontsize = 8, col = 'black', lwd = 10),ch = 0, axLab.pos = 0.06, rot = c(0, -60, 60, c(rep(0, nrow(fread(nds))))))
  dev.off()
  bla <- fread(nds)
  blax <- merge(func_ann$names[.(bla$node.lab),,nomatch=0], bla[,list(node.lab, radius)], by.x = 'ID', by.y = 'node.lab')[order(radius),]
  blay <- merge(bla, hpd$nodes, by.x = 'node.lab', by.y = 'lab')
  setkey(blay, axis, radius.y)
  blay[,term:=func_ann$names[.(node.lab),Name]]
  blay[,list(term,node.lab,node.text,radius.y,radius.x)]
  return(blay[,list(term,node.lab,node.text,radius.y,radius.x)])
})

blub <- rbindlist(lapply(cols[1:9], function(i) {
  cat(i, '\n')
  tmp <- hpd
  # tmp$edges[!tmp$edges$color%in%i,'color'] <- bgcolor # color all edges not of color i in the background color
  # tmp$edges <- tmp$edges[c(which(!tmp$edges$color%in%bgcolor)),]
  tmp$edges <- tmp$edges[tmp$edges$color%in%i,] # select edges of color i
  if (nrow(tmp$edges)==0)
    return(data.table(lab1=character(0),Name=character(0),lab2=character(0),Tissue=character(0)))
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
  require(colorspace)
  labs <- names(.cs$tooColor)
  names(labs) <- col2hcl(.cs$tooColor)
  res[,Tissue:=labs[hpd$nodes[hpd$nodes$axis%in%3,"lab"][hpd$nodes[hpd$nodes$axis%in%3,'color']==i]]]
  return(res)
}))

setkey(blub, lab2)
setkey(modulemembers, Module)

res <- blub[modulemembers,,allow.cartesian=T,nomatch=0]
setnames(res, c('Gene set ID', 'Gene set name', 'Module', 'Tissue of origin', 'Gene name/p-site'))
# setnames(res, c('ID', 'Name', 'WGCNA Cluster', 'Protein/p-site'))
fwrite(res, '~/phospho/Chen/Res/20180523_moreSuppTables/WGCNA_new_annotations_NCI60_new.txt', col.names = T, row.names = F, sep = '\t')
