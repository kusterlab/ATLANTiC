source('~/atlantic/shiny/ilandscape/landscape_functions.R')
source('~/phospho/Chen/pubFigures2/R/colorScheme.R')

corandp <- function (x, y = NULL, use = "pairwise.complete.obs", alternative = c("two.sided", "less", "greater"), nThreads = 16, ...) 
{
  ia = match.arg(alternative)
  cat('Correlation\n')
  cor = WGCNA:::cor(x, y, use = use, nThreads = nThreads, ...)
  x = as.matrix(x)
  finMat = !is.na(x)
  if (is.null(y)) {
    cat('Observations\n')
    np = t(finMat) %*% finMat
  }
  else {
    cat('Observations\n')
    y = as.matrix(y)
    np = t(finMat) %*% (!is.na(y))
  }
  if (ia == "two.sided") {
    cat('P-Value\n')
    T = sqrt(np - 2) * abs(cor)/sqrt(1 - cor^2)
    # cor <- lowerTriangle(cor)
    gc()
    p = 2 * pt(T, np - 2, lower.tail = FALSE)
  }
  else if (ia == "less") {
    cat('P-Value\n')
    T = sqrt(np - 2) * cor/sqrt(1 - cor^2)
    # cor <- lowerTriangle(cor)
    gc()
    p = pt(T, np - 2, lower.tail = TRUE)
  }
  else if (ia == "greater") {
    cat('P-Value\n')
    T = sqrt(np - 2) * cor/sqrt(1 - cor^2)
    # cor <- lowerTriangle(cor)
    gc()
    p = pt(T, np - 2, lower.tail = FALSE)
  }
  list(cor = cor, p = p, nObs = np)
}

# # WGCNA NCI60
# 
# mat <- readRDS(file = "/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180503_gsva/gsva_nci60_mods.RDS")
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/landscape_map_nci60_wgcna.png', width = (210-35*2), height = (210-35*2), units = 'mm', pointsize = 10, res = 300)
# par(mar=c(6.1, 3.1, 2.1, 1.1))
# landscape(mat = mat, cutfac = 1, sealevel = 0.8, sub = '_NCI60', main = 'Activity landscape of the NCI60 dataset - wgcna')
# dev.off()

# PATH NCI60

nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist.RDS')
for (i in seq_along(nci60)) {
  colnames(nci60[[i]]) <- gsub('_NCI60', '', colnames(nci60[[i]]))
}
names(nci60) <- c('NetworKIN', 'Pathways', 'Integrated', 'SIGNOR', 'Integrated_NetworKIN')
# 
# termdis <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/msigdb.rds')
# termdis <- as.dist(as.matrix(termdis)[rownames(nci60$path),rownames(nci60$path)])
# require(dendsort)
# term_clust <- dendsort(hclust(termdis, method = 'ward.D2'), type = 'min')
# term_order <- term_clust$labels[term_clust$order]

nci60_path <- clu(nci60$Pathways, dist.method.row = 'correlation', dist.method.col = 'euclidean', clust.method.row =  'mcquitty', clust.method.col =  'mcquitty')
colnames(nci60_path) <- gsub('_NCI60', '', colnames(nci60_path))
# nci60_path <- nci60_path[term_order,]
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/landscape_map_nci60_path.png', width = (210-35*2), height = (210-35*2), units = 'mm', pointsize = 10, res = 300)
# par(mar=c(6.1, 3.1, 2.1, 1.1))
# landscape(mat = nci60$path, cutfac = 1, sealevel = 0.8, sub = '_NCI60', main = 'Activity landscape of the NCI60 dataset - pathways')
# dev.off()

# # NKIN NCI60
# 
# nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist.RDS')
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/landscape_map_nci60_nkin.png', width = (210-35*2), height = (210-35*2), units = 'mm', pointsize = 10, res = 300)
# par(mar=c(6.1, 3.1, 2.1, 1.1))
# landscape(mat = nci60$nkin, cutfac = 1, sealevel = 0.8, sub = '_NCI60', main = 'Activity landscape of the NCI60 dataset - nkin')
# dev.off()

# # INTEG NCI60
# 
# nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist.RDS')
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/landscape_map_nci60_integ.png', width = (210-35*2), height = (210-35*2), units = 'mm', pointsize = 10, res = 300)
# par(mar=c(6.1, 3.1, 2.1, 1.1))
# landscape(mat = nci60$integ, cutfac = 1, sealevel = 0.8, sub = '_NCI60', main = 'Activity landscape of the NCI60 dataset - integ')
# dev.off()

# # SIGNOR NCI60
# 
# nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist.RDS')
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/landscape_map_nci60_signor.png', width = (210-35*2), height = (210-35*2), units = 'mm', pointsize = 10, res = 300)
# par(mar=c(6.1, 3.1, 2.1, 1.1))
# landscape(mat = nci60$signor, cutfac = 1, sealevel = 0.8, sub = '_NCI60', main = 'Activity landscape of the NCI60 dataset - signor')
# dev.off()



# # WGCNA CRC65
# 
# mat <- readRDS(file = "/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180503_gsva/gsva_crc65_mods.RDS")
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/landscape_map_crc65_wgcna.png', width = (210-35*2), height = (210-35*2), units = 'mm', pointsize = 10, res = 300)
# par(mar=c(6.1, 3.1, 2.1, 1.1))
# res <- landscape(mat = mat, cutfac = 1, sealevel = 0.8, sub = '_CRC65', main = 'Activity landscape of the CRC65 dataset - wgcna')
# dev.off()

# PATH CRC65

crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist_crc65.RDS')
for (i in seq_along(crc65)) {
  colnames(crc65[[i]]) <- gsub('_CRC65', '', colnames(crc65[[i]]))
}
names(crc65) <- c('NetworKIN', 'Pathways', 'Integrated', 'SIGNOR', 'Integrated_NetworKIN')

crc65_path <- clu(crc65$Pathways, dist.method.row = 'correlation', dist.method.col = 'euclidean', clust.method.row =  'mcquitty', clust.method.col =  'mcquitty')
colnames(crc65_path) <- gsub('_CRC65', '', colnames(crc65_path))
# crc65_path <- crc65_path[term_order,]

path <- list(CRC65=crc65, NCI60=nci60)

path <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/zscore_data_new.RDS')
saveRDS(path, '/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/clustered_data.RDS', compress = T)

selection <- list(
  list(col=c(57:64),row=c(589:623)),
  list(col=c(54:61),row=c(748:770)),
  list(col=c(48:55),row=c(811:845)),
  list(col=c(22:30),row=c(994:1025)),
  list(col=c(11:14),row=c(935:984)),
  list(col=c(1:8),row=c(652:737)),
  list(col=c(23:28),row=c(91:128))
)
setEPS()
postscript("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigure3/figures/landscape_crc65_new.eps", width = 8.7/2.539998, height = 5.8/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular')
res_crc65 <- landscape(mat = path$CRC65$Pathways, cutfac = 1, sealevel = 0.35, selection = selection, main = 'Activity landscape of the CRC65 dataset', bw = 2.5, dx = 0.5, dy = 1.5, cex.axis = 1, mar = c(1.5, 1.5, 1.5, 1.1), ylab = 'Pathways')
dev.off()

selection <- list(
  list(col=c(53:60),row=c(1055:1100)),
  list(col=c(19:23),row=c(1106:1133)),
  list(col=c(46:50),row=c(32:73)),
  list(col=c(33:40),row=c(664:709)),
  list(col=c(53:60),row=c(1174:1233)),
  list(col=c(34:38),row=c(209:248))
)

setEPS()
postscript("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigure3/figures/landscape_nci60_path_new.eps", width = 10.5/2.539998, height = 7/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular')
res_nci60_path <- landscape(mat = path$NCI60$Pathways, cutfac = 1, sealevel = 0.35, selection = selection, main = 'Activity landscape of the NCI60 dataset', bw = 2.5, dx = 0.5, dy = 1.5, cex.axis = 1, mar = c(1.5, 1.5, 1.5, 1.1), ylab = 'Pathways')
dev.off()


matchobject <- colnames(res_crc65$mat)
source('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/src/crc65_doublingtime.R')

overl <- intersect(matchobject, rownames(FUannotation))

bla <- corandp(t(res_crc65$mat[unique(res_crc65$land[region=='F',pathway]),overl]),FUannotation[overl,"Lag time (days)",drop=F], use = 'p')
FUannotation[res_crc65$land[region=='F',unique(cellline)],"Doubling time (hours)",drop=F]

cor.test(colMeans(res_crc65$mat[grep('CELL_CYCLE', unique(res_crc65$land[region=='F',pathway]), ignore.case = T, value = T),overl], na.rm = T),FUannotation[overl,"Doubling time (hours)"], use = 'p')

selection <- list(list(col=c(3:3),row=c(1:235)),
                  list(col=c(4:4),row=c(1:235)),
                  list(col=c(11:17),row=c(115:134)),
                  list(col=c(29:34),row=c(30:51)),
                  list(col=c(47:54),row=c(66:80)),
                  list(col=c(58:60),row=c(90:111))
)

setEPS()
postscript("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigure3/figures/landscape_nci60_new.eps", width = 8.7/2.539998, height = 5.8/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular')
res_nci60 <- landscape(mat = path$NCI60$`Kinases (by ATLANTiC Score)`, cutfac = 1, sealevel = 0.35, selection = selection, main = 'Activity landscape of the NCI60 dataset', bw = 2.5, dx = 0.5, dy = 1.5, cex.axis = 1, mar = c(1.5, 1.5, 1.5, 1.1), ylab = 'Kinases (by ATLANTIC Score)')
dev.off()

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

# plot(sort(res_nci60$smat[,"K562"], decreasing = T), axes = F, xlab = '', ylab = '', pch = 19, cex = 0.7, col = map2color(sort(res_nci60$smat[,"K562"], decreasing = T), pal = res_nci60$col), ylim = c(0,1))
# par(new = T)
postscript("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigure3/figures/waterfall_K562.eps", width = 8.7/2.539998, height = 5.8/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular')
par(mar=c(2.1, 4.1, 2.1, 2.1))
plot(sort(res_nci60$mat[,"K562"], decreasing = T), axes = F, xlab = '', ylab = '', pch = 19, cex = 0.7, col = map2color(sort(res_nci60$smat[,"K562"], decreasing = T), pal = res_nci60$col), ylim = c(0,1))
box()
# axis(1)
axis(2, las = 2)
text(x = c(1,2,3)+20,y = sort(res_nci60$mat[,"K562"], decreasing = T)[c(1,2,3)], labels = names(sort(res_nci60$mat[,"K562"], decreasing = T)[c(1,2,3)]))
mtext(text = 'Ranked kinase activity - K562', side = 3, line = 0.5)
mtext(text = 'Relative activity', side = 2, line = 3)
mtext(text = 'Kinases', side = 1, line = 0.5)
dev.off()

postscript("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigure3/figures/waterfall_SR.eps", width = 8.7/2.539998, height = 5.8/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular')
par(mar=c(2.1, 4.1, 2.1, 2.1))
plot(sort(res_nci60$mat[,"SR"], decreasing = T), axes = F, xlab = '', ylab = '', pch = 19, cex = 0.7, col = map2color(sort(res_nci60$smat[,"SR"], decreasing = T), pal = res_nci60$col), ylim = c(0,1))
box()
# axis(1)
axis(2, las = 2)
text(x = c(1,2,3)+20,y = sort(res_nci60$mat[,"SR"], decreasing = T)[c(1,2,3)], labels = names(sort(res_nci60$mat[,"SR"], decreasing = T)[c(1,2,3)]))
mtext(text = 'SR', side = 3, line = 0.5)
mtext(text = 'Relative activity', side = 2, line = 3)
mtext(text = 'Kinases', side = 1, line = 0.5)
dev.off()


selection <- list(
  list(col=c(26:26),row=c(1:256)),
  list(col=c(15:15),row=c(1:256)),
  list(col=c(17:17),row=c(1:256)),
  list(col=c(61:64),row=c(144:165)),
  list(col=c(46:55),row=c(71:91))
)

setEPS()
postscript("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigure3/figures/landscape_crc65_atlantic_new.eps", width = 10.5/2.539998, height = 7/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular')
res_crc65_atlantic <- landscape(mat = path$CRC65$`Kinases (by ATLANTiC Score)`, cutfac = 1, sealevel = 0.35, selection = selection, main = 'Activity landscape of the CRC65 dataset', bw = 2.5, dx = 0.5, dy = 1.5, cex.axis = 1, mar = c(1.5, 1.5, 1.5, 1.1), ylab = 'Kinases (by ATLANTIC Score)')
dev.off()

plotPath <- function(land, area, topN = NULL, file = './example.png') {
  if(is.null(topN)) {topN <- nrow(land)}
  tmp <- land[region==area,list(`Relative activity`=sum(score)),by=pathway][order(`Relative activity`, decreasing = T),][1:topN,][!is.na(pathway),][,`Relative activity`:=`Relative activity`/max(`Relative activity`)]
  setkey(tmp, `Relative activity`)
  labs <- gsub('c\\("|"|\\)|,', '', str_wrap(str_split(tmp[,pathway], '_'), width = 40))
  nr <- sum(str_count(labs, '\n')+1)
  png(filename = file, width = (210-35*2), height = nr*10, units = 'mm', pointsize = 10, res = 300)
  par(mar=c(2.1,20,3.1,1.1))
  tmp[,barplot(`Relative activity`, names.arg = labs, horiz = T, las = 2, cex.names = 1, col = .cs$tumblue, axes = F, yaxs = 'i')]  
  axis(side = 3, xaxs = 'i', cex.axis = 1, yaxs = 'i')
  mtext(side = 3, text = 'Relative importance', line = 2, cex = 1)
  dev.off()
}

file <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/A_test.png'
plotPath(land = res$land, area = 'A', topN = NULL, file = file)
file <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/B_test.png'
plotPath(land = res$land, area = 'B', topN = NULL, file = file)
file <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/C_test.png'
plotPath(land = res$land, area = 'C', topN = NULL, file = file)
file <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/D_test.png'
plotPath(land = res$land, area = 'D', topN = NULL, file = file)
file <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/E_test.png'
plotPath(land = res$land, area = 'E', topN = NULL, file = file)
file <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/F_test.png'
plotPath(land = res$land, area = 'F', topN = NULL, file = file)
file <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/G_test.png'
plotPath(land = res$land, area = 'G', topN = NULL, file = file)




res$land[region=='A',sum(score),by=list(cellline,pathway)][cellline==unique(cellline)[4],barplot(sort(V1), names.arg = pathway, horiz = T, las = 2, cex.names = 0.6, main = unique(cellline))]
res$land[region=='B',list(sum(smoothscore),sum(score)),by=pathway][,cor(V1,V2)]

# # NKIN CRC65
# 
# crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist_crc65.RDS')
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/landscape_map_crc65_nkin.png', width = (210-35*2), height = (210-35*2), units = 'mm', pointsize = 10, res = 300)
# par(mar=c(6.1, 3.1, 2.1, 1.1))
# landscape(mat = crc65$nkin, cutfac = 1, sealevel = 0.8, sub = '_CRC65', main = 'Activity landscape of the CRC65 dataset - nkin')
# dev.off()
# 
# # INTEG CRC65
# 
# crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist_crc65.RDS')
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/landscape_map_crc65_integ.png', width = (210-35*2), height = (210-35*2), units = 'mm', pointsize = 10, res = 300)
# par(mar=c(6.1, 3.1, 2.1, 1.1))
# landscape(mat = crc65$integ, cutfac = 1, sealevel = 0.8, sub = '_CRC65', main = 'Activity landscape of the CRC65 dataset - integ')
# dev.off()
# 
# # SIGNOR CRC65
# 
# crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/Res/20180625_kinaseActivity/scorelist_crc65.RDS')
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/landscape_map_crc65_signor.png', width = (210-35*2), height = (210-35*2), units = 'mm', pointsize = 10, res = 300)
# par(mar=c(6.1, 3.1, 2.1, 1.1))
# landscape(mat = crc65$signor, cutfac = 1, sealevel = 0.8, sub = '_CRC65', main = 'Activity landscape of the CRC65 dataset - signor')
# dev.off()

















# pheatmap(smt, cluster_rows = F, cluster_cols = F, color = col, gaps_row = seq(20,nrow(bla),by = 20), gaps_col = seq(6, ncol(bla), by = 6))
# pth <- image.smooth(smt, theta = 3)$x
# cll <- image.smooth(smt, theta = 3)$y
# act <- image.smooth(smt, theta = 3)$z
# png(filename = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/landscape.png', width = (210-35*2), height = (210-35*2), units = 'mm', pointsize = 10, res = 300)
# par(mar = c(0,2,0,0))
# persp(x = pth, y = cll, z = act, theta = 45, phi = 30, r = 10000, d = 10, box = T, border = 'grey30', shade = NA, xlab = 'Pathways', ylab = 'Celllines', zlab = 'Activity')
# dev.off()
# z <- c(smt)
# bla <- smt
# colnames(bla) <- x
# rownames(bla) <- y
# bla <- melt(bla)
# colnames(bla) <- c('x','y','z')
# scatterplot3d(x = bla$x, y = bla$y, z = bla$z, pch = '.')
# bla <- list(x = x, y = y, z = z)
# persp(bla, phi = 30, theta = 20, d = 5)
# melt(smt[rord,cord])
# bla <- kde2d(smt)
# sm.density(melt(smt), h=c(5,5), nbins=361,
#            xlim=c(-180,180), ylim=c(-180,180),
#            xlab="phi", ylab="psi", zlab="")
