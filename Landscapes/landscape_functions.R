require(shiny)
require(sm)
require(MASS)
require(scales)
require(fields)
require(plyr)
require(data.table)
require(stringr)
require(parallelDist)
require(WGCNA)

geoMean <- function (x, na.rm = TRUE) 
{
  if (is.null(nrow(x))) {
    exp(mean(log(x), na.rm = TRUE))
  }
  else {
    exp(apply(log(x), 2, mean, na.rm = na.rm))
  }
}
geoSD <- function (x, na.rm = TRUE) 
{
  if (is.null(nrow(x))) {
    exp(sd(log(x), na.rm = T))
  }
  else {
    exp(apply(log(x), 2, sd, na.rm = na.rm))
  }
}

clu <- function(mat,
                dist.method.row = c('correlation',
                                    'euclidean',
                                    'maximum',
                                    'manhattan',
                                    'canberra',
                                    'binary',
                                    'minkowski'),
                dist.method.col = c('correlation',
                                    'euclidean',
                                    'maximum',
                                    'manhattan',
                                    'canberra',
                                    'binary',
                                    'minkowski'),
                clust.method.row = c('ward.D',
                                     'ward.D2',
                                     'single',
                                     'complete',
                                     'average',
                                     'mcquitty',
                                     'median',
                                     'centroid'),
                clust.method.col = c('ward.D',
                                     'ward.D2',
                                     'single',
                                     'complete',
                                     'average',
                                     'mcquitty',
                                     'median',
                                     'centroid')) {
  dist.method.row <- match.arg(dist.method.row)
  dist.method.col <- match.arg(dist.method.col)
  clust.method.row <- match.arg(clust.method.row)
  clust.method.col <- match.arg(clust.method.col)
  # mat <- t(scale(t(mat)))
  # mat <- t(apply(mat, 1, function(i) (i-geoMean(i, na.r = T))/geoSD(i, na.rm = T)))
  
  mat <- rescale(mat, to = c(0,1))
  
  # smt <- t(apply(mat, 1, rescale, c(0,1)))
  smt <- mat
  # Rows
  if (dist.method.row=='correlation') {
    rord <- fastcluster:::hclust(as.dist(1-cor(t(smt), use = 'p', method = 'pearson', nThreads = 10, quick = 1)), method = clust.method.row)
    rord$labels <- rownames(smt)
  } else {
    rord <- fastcluster:::hclust(parDist(smt, method = dist.method.row, threads = 10), method = clust.method.row)
    rord$labels <- rownames(smt)
  }
  # Cols
  if (dist.method.col=='correlation') {
    cord <- fastcluster:::hclust(as.dist(1-cor(smt, use = 'p', method = 'pearson', nThreads = 10, quick = 1)), method = clust.method.col)
    cord$labels <- colnames(smt)
  } else {
    cord <- fastcluster:::hclust(parDist(t(smt), method = dist.method.col, threads = 10), method = clust.method.col)
    cord$labels <- colnames(smt)
  }
  cord <- cord$labels[cord$order]
  rord <- rord$labels[rord$order]
  smt <- smt[rord, cord] #
  return(smt)
}

smo <- function(smt, cutfac = 1, bw = 3, dx = 1, dy = 1) {
  if(bw==0) {
    bw <- 1e-13
  }
  if(dx==0) {
    dx <- 1e-13
  }
  if(dy==0) {
    dy <- 1e-13
  }
  cut <- sd(smt)*cutfac
  med <- median(smt)
  smt[smt>(med+cut)] <- med+cut
  smt[smt<(med-cut)] <- med-cut
  bla <- image.smooth(smt, theta = bw, xwidth = NULL, ywidth = NULL, dx = dx, dy = dy)$z
  rownames(bla) <- rownames(smt)
  colnames(bla) <- colnames(smt)
  bla <- rescale(bla, c(0,1))
  return(bla)
}

sea <- function(bla, sealevel = 0.8, ncol = 256) {
  t1 <- table(seq(min(bla),max(bla), length.out = ncol)>=quantile(bla, sealevel))
  landfrac <- t1["TRUE"]/sum(t1)
  col <- colorRampPalette(c(colorRampPalette(rev(blues9[3:9]))(ncol*(1-landfrac)/10), colorRampPalette(terrain.colors(n = 10))(ncol*landfrac/10)))(ncol)
  return(col)
}

getLand <- function(smt, bla, hil, label = 'pathway', sealevel = 0.8) {
  smtdt <- setDT(melt(smt))
  setnames(smtdt, c(label, 'cellline', 'score'))
  setkeyv(smtdt, c(label, 'cellline'))
  bladt <- setDT(melt(bla))
  setnames(bladt, c(label, 'cellline', 'smoothscore'))
  setkeyv(bladt, c(label, 'cellline'))
  land <- data.table(which(bla>quantile(bla, sealevel), arr.ind = T), keep.rownames = T)
  setnames(land, 'rn', label)
  land[,cellline:=colnames(bla)[col]]
  setkeyv(land, c(label, 'cellline'))
  land <- bladt[land]
  setkeyv(land, c(label, 'cellline'))
  land <- smtdt[land]
  setkeyv(land, c(label, 'cellline'))
  pltdt <- setDT(melt(hil))
  setnames(pltdt, c('cellline', label, 'region'))
  pltdt[,region:=LETTERS[region-1]]
  setkeyv(pltdt, c(label, 'cellline'))
  land <- pltdt[land]
  land[,eval(label):=as.character(eval(as.name(label)))]
  land[,cellline:=as.character(cellline)]
  land[,score:=round(score, digits = 3)]
  land[,smoothscore:=round(smoothscore, digits = 3)]
  setkey(land, row, col, score)
  pthw <- land[!is.na(region),][,list(`Relative activity`=sum(score)),by=list(region, eval(as.name(label)))][order(`Relative activity`, decreasing = T),]
  setnames(pthw, 'as.name', label)
  pthw[,`Relative activity`:=`Relative activity`/max(`Relative activity`),by=list(region)]
  pthw[,`Relative activity`:=round(`Relative activity`, 2)]
  # setnames(pthw, 'pathway', 'Pathway')
  res <- list(land = land, pathway = pthw)
	names(res) <- c('land', label)
  return(res)
}

reorderMat <- function(bla) {
  plt <- t(bla[rev(rownames(bla)),])
  return(plt)
}

getRegion <- function(plt, coords) {
  pc <- coords[c('xmin', 'xmax', 'ymin', 'ymax')]
  if(!is.null(pc)) {
    pc <- lapply(pc, as.numeric)
    pc <- lapply(pc, round, 0)
    region <- list(list(col=c(pc$xmin:pc$xmax),row=c((ncol(plt)-pc$ymin+1):(ncol(plt)-pc$ymax+1))))
    region[[1]] <- lapply(region[[1]], function(i) i[i>0])
    region[[1]]$row <-  region[[1]]$row[region[[1]]$row<=ncol(plt)]
    region[[1]]$col <-  region[[1]]$col[region[[1]]$col<=nrow(plt)]
  } else {
    region <- list()
  }
  return(region)
}

getRow <- function(plt, coords) {
  # column() <- NULL
  pc <- coords[c('y')]
  if(!is.null(pc)) {
    pc <- lapply(pc, as.numeric)
    pc <- lapply(pc, round, 0)
    ROW <- list(list(col=c(1:nrow(plt)),row=c(ncol(plt)-pc$y+1)))
    ROW[[1]] <- lapply(ROW[[1]], function(i) i[i>0])
    ROW[[1]]$row <-  ROW[[1]]$row[ROW[[1]]$row<=ncol(plt)]
    ROW[[1]]$col <-  ROW[[1]]$col[ROW[[1]]$col<=nrow(plt)]
  } else {
    ROW <- list()
  }
  return(ROW)
}

getColumn <- function(plt, coords) {
  # row() <- NULL
  pc <- coords[c('x')]
  if(!is.null(pc)) {
    pc <- lapply(pc, as.numeric)
    pc <- lapply(pc, round, 0)
    COL <- list(list(col=c(pc$x),row=c(1:ncol(plt))))
    COL[[1]] <- lapply(COL[[1]], function(i) i[i>0])
    COL[[1]]$row <-  COL[[1]]$row[COL[[1]]$row<=ncol(plt)]
    COL[[1]]$col <-  COL[[1]]$col[COL[[1]]$col<=nrow(plt)]
  } else {
    COL <- list()
  }
  return(COL)
}

decideSelection <- function(row, column, region) {
  if (!is.null(region)) {
    return(region)
  }
  if (!is.null(column)) {
      return(column)
  }
  if (!is.null(row)) {
    return(row)
  }
  return(list())
}

restrictPlt <- function(plt, sealevel = 0.8, selection = NULL) {
  if (length(selection)>26) {
    warning('More than 26 regions selected! Only the first 26 regions will be shown!')
    selection <- selection[1:26]
  }
  reg <- selection
  cut <- quantile(plt, sealevel)
  rct <- list()
  for (i in seq_along(reg)) {
    i <- as.integer(i)
    plt[reg[[i]]$col,rev(colnames(plt))[reg[[i]]$row]][plt[reg[[i]]$col,rev(colnames(plt))[reg[[i]]$row]]>=cut] <- i+1L
    rct[[i]] <- list(xleft=min(reg[[i]]$col)-0.5,
                     ybottom=match(rev(colnames(plt))[max(reg[[i]]$row)],colnames(plt))-0.5,
                     xright=max(reg[[i]]$col)+0.5,
                     ytop=match(rev(colnames(plt))[min(reg[[i]]$row)],colnames(plt))+0.5)
  }
  plt[plt<2L] <- 28L
  hil <- plt
  res <- list(hil=hil, rct=rct)
  return(res)
}

plotMap <- function(plt, res, col, main, cex.axis = 0.6, cex.lab = 1, mar = c(7.1, 4, 2.5, 0.6), ylab = 'Pathways') {
  par(mar=mar)
  image.plot(plt, col = col, x = 1:nrow(plt), y = 1:ncol(plt), axes = F, xlab = NA, ylab = NA, legend.shrink = 1, legend.lab = 'Relative activity', legend.width = 1, legend.line = 3, legend.mar = 6.1, legend.cex = cex.lab, axis.args = list(cex.axis = cex.axis, tck = -0.2, mgp = c(3,0.3,0)))
  box()
  abline(v = seq(0,round_any(nrow(plt), accuracy = 10), length.out = 6)+0)
  abline(h = seq(0,round_any(ncol(plt), accuracy = 10), length.out = 6)+0)
  Map(function(x,y,z) 
    axis(1,at=x,col.axis=y,labels=z,lwd=0,las=2, cex.axis = cex.axis, mgp = c(3,0.3,0)),
    1:nrow(plt),
    c(1,1)[rownames(plt)%in%c('NCI-H716', 'C10', 'K562', 'SR', 'SW620', 'SW480', 'LS 174T', 'LS 180', 'HDC-57', 'HDC-54')+1],
    rownames(plt)
  )
  axis(1, at = 1:nrow(plt), labels = F, tck = -0.005)
  axis(2, at = ceiling(seq(round_any(ncol(plt), accuracy = 10), 1, length.out = 20)), labels = rev(ceiling(seq(round_any(ncol(plt), accuracy = 10), 1, length.out = 20))), las = 2, cex.axis = cex.axis, tck = -0.005, mgp = c(3,0.3,0))
  mtext(text = 'Cell lines', side = 1, line = 6, cex = cex.lab)
  mtext(text = ylab, side = 2, line = 2.5, cex = cex.lab)
  mtext(text = main, side = 3, line = 1, cex = cex.lab)
  hicol <- alpha('white', 0)
  dimcol <- alpha('black', 0.3)
  cols <- c(rep(hicol, times = length(unique(c(res$hil)))-1), dimcol)
  if (length(cols) == 1) cols <- append(hicol, cols)
  if(!all(res$hil==28)) {
    image.plot(res$hil, col = cols, x = 1:nrow(res$hil), y = 1:ncol(res$hil), add = T, legend.mar = 4.1, legend.shrink = 1, legend.width = 0.0000001, legend.lab = NA, axis.args = list(line = 1e6))
    for (i in seq_along(res$rct)) {
      i <- as.integer(i)
      rect(xleft = res$rct[[i]]$xleft, ybottom = res$rct[[i]]$ybottom, xright = res$rct[[i]]$xright, ytop = res$rct[[i]]$ytop, border = 'black', lty = 2, lwd = 1.5)
      if (length(res$rct)>1) {
        text(x = res$rct[[i]]$xleft+((res$rct[[i]]$xright)-(res$rct[[i]]$xleft))/2, y = res$rct[[i]]$ybottom+((res$rct[[i]]$ytop)-(res$rct[[i]]$ybottom))/2, labels = LETTERS[i], adj = c(0.5,0.5), col = 'black', font = 2)
      }
    }
  }
}

plotMap2 <- function(plt, res, col, main, cex.axis = 0.6, cex.lab = 1, mar = c(2.1, 2.1, 2.1, 1.1), ylab = 'Pathways') {
  par(mar=mar)
  image.plot(plt, col = col, x = 1:nrow(plt), y = 1:ncol(plt), axes = F, xlab = NA, ylab = NA, legend.shrink = 1, legend.lab = 'Relative activity', legend.width = 1, legend.line = 2, legend.mar = 4.1, legend.cex = cex.lab, axis.args = list(cex.axis = cex.axis, tck = -0.2, mgp = c(3,0.3,0)))
  box()
  abline(v = seq(0,round_any(nrow(plt), accuracy = 10), length.out = 6)+0.5, lwd = 0)
  abline(h = seq(0,round_any(ncol(plt), accuracy = 10), length.out = 6)+0.5, lwd = 0)
  mtext(text = 'Cell lines', side = 1, line = 0.5, cex = cex.lab)
  mtext(text = ylab, side = 2, line = 0.5, cex = cex.lab)
  mtext(text = main, side = 3, line = 0.5, cex = cex.lab)
  hicol <- alpha('white', 0)
  dimcol <- alpha('black', 0.3)
  cols <- c(rep(hicol, times = length(unique(c(res$hil)))-1), dimcol)
  if (length(cols) == 1) cols <- append(hicol, cols)
  if(!all(res$hil==28)) {
    for (i in seq_along(res$rct)) {
      i <- as.integer(i)
      rect(xleft = res$rct[[i]]$xleft, ybottom = res$rct[[i]]$ybottom, xright = res$rct[[i]]$xright, ytop = res$rct[[i]]$ytop, border = 1, lty = 1, lwd = 0.75)
      text(x = res$rct[[i]]$xleft+((res$rct[[i]]$xright)-(res$rct[[i]]$xleft))/2, y = res$rct[[i]]$ybottom+((res$rct[[i]]$ytop)-(res$rct[[i]]$ybottom))/2, labels = LETTERS[i], adj = c(0.5,0.5), col = 1, font = 2)
    }
  }
}

landscape <- function(mat,
                      cutfac = 1,
                      bw = 3,
                      dx = 1,
                      dy = 1,
                      sealevel = 0.8,
                      dist.method.row = 'correlation',
                      dist.method.col = 'euclidean',
                      clust.method.row = 'mcquitty',
                      clust.method.col = 'mcquitty',
                      selection = NULL, 
                      main,
                      ylab,
                      cex.axis = 0.6,
                      mar = c(6.1, 3.5, 2.5, 1.1)) {
  smt <- clu(mat, dist.method.row = dist.method.row, dist.method.col = dist.method.col, clust.method.row = clust.method.row, clust.method.col = clust.method.col)
  bla <- smo(smt, cutfac = cutfac, bw = bw, dx = dx, dy = dy)
  col <- sea(bla, sealevel = sealevel)
  plt <- reorderMat(bla)
  res <- restrictPlt(plt, sealevel = sealevel, selection = selection)
  plotMap2(plt, res, col = col, main = main, ylab = ylab, cex.axis = cex.axis, mar = mar)
  res <- getLand(smt = smt, bla = bla, hil = res$hil, sealevel = sealevel)
  return(list(mat = smt, smat = bla, land = res$land, path = res$pathway, col = col))
}
