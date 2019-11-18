require(reshape2)
source('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/src/20170809_fisher.R')
source('/media/kusterlab/internal_projects/active/CRC65/data/24_CRC64_1and2and3_120ppm_1.4.1.2_peakproperties/analysis/colors.R')
source('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/src/multiSetMEs_FIX.R')

compareGroupings <- function(wgcna1, wgcna2) {
  wgcna2[,present:=T]
  setkey(wgcna1, id)
  anntab <- as.data.frame(acast(data = wgcna2, formula = id~moduleColor, value.var = 'present', fill = F))
  anntab$moduleColor <- wgcna1[.(rownames(anntab)),moduleColor] # HERE!
  fisher.rbind(subtypename = 'moduleColor',
               categories = colnames(anntab)[!colnames(anntab)%in%'moduleColor'],
               func_ann = NULL,
               data = anntab,
               alternative = 'two.sided',
               method.subtype = 'BH',
               method.all = 'none',
               digits = 3,
               cutoff = 1,
               multimembersep = NULL,
               subtypestoexclude = 'background',
               go_rollup = NULL,
               ncores = 10,
               parseResult = T)
}

alldup <- function (value) { 
  duplicated(value) | duplicated(value, fromLast = TRUE)
}

nobsfun <- function(x) {
  finMat = !is.na(x)
  nobs <- t(finMat) %*% finMat
  return(nobs)
}

generateTermAssignment <- function(ids, func_ann, go_rollup=func_ann$go_offspring, identifier = 'SYMBOL') {
  require(pbmcapply)
  setkeyv(func_ann$annotations, identifier)
  nongo <- func_ann$annotations[.(ids),,nomatch = 0][!grepl('GO_', ONTOLOGY),list(eval(as.name(identifier)), ID)]
  setnames(nongo, 'V1', identifier)
  setkey(nongo)
  nongo <- unique(nongo)
  nongo[,PID:=ID]
  godt <- data.table(Parent=rep(names(go_rollup), sapply(go_rollup, length)),
                     Child=unname(unlist(go_rollup)))
  godt <- rbindlist(list(godt, data.table(Parent=godt[!is.na(Child),unique(Parent)],Child=godt[!is.na(Child),unique(Parent)])))
  godt[is.na(Child),Child:=Parent]
  go <- func_ann$annotations[.(ids),,nomatch = 0][grepl('GO_', ONTOLOGY),list(eval(as.name(identifier)), ID)]
  setnames(go, 'V1', identifier)
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

plotManyTraces <- function(mat, idlist, sampleid='Samples', featureid='Proteins', valueid='Expression', plotsperpage = 10, file = 'Traces.pdf') {
  require(gridExtra)
  plotlist <- lapply(seq_along(idlist), function(i) plotTraces(mat = mat, ids = idlist[[i]], sampleid = sampleid, featureid = featureid, valueid = valueid, title = names(idlist)[i]))
  names(plotlist) <- names(idlist)
  ncol <- 2
  pl <- marrangeGrob(grobs = plotlist, nrow=plotsperpage/ncol, ncol=ncol)
  ggsave(file,
         pl,
         device = 'pdf', 
         width = 210, 
         height = 297, 
         units = "mm", onefile = T)
}

plotTraces <- function(mat, ids, sampleid='Samples', featureid='Proteins', valueid='Expression', title = 'Traces') {
  require(reshape2)
  require(ggplot2)
  matdt <- data.table(melt(mat))
  setnames(matdt, c(featureid, sampleid, valueid))
  if (is.factor(matdt[,eval(as.name(sampleid))])) {
    matdt[,c(sampleid):=factor(eval(as.name(sampleid)), levels = colnames(mat), labels = colnames(mat), ordered = T)]
  }
  # matdt[eval(as.name(featureid))%in%ids,highlight:='high']
  # matdt[is.na(highlight),highlight:='low']
  pl <- ggplot(data = matdt[eval(as.name(featureid))%in%ids,], aes(x = eval(as.name(sampleid)), y = eval(as.name(valueid)), colour = eval(as.name(featureid)), group = eval(as.name(featureid)))) + # , alpha = highlight
    geom_line() +
    theme_classic() +
    theme(legend.position = 'right',
          axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(size = rel(0.7)),
          axis.title = element_text(size = rel(0.7))) +
    labs(title = title, x = sampleid, y = valueid) +
    scale_color_discrete(name = featureid) +
    # scale_y_continuous(labels=trans_format('log10',math_format(10^.x))) +
    # scale_y_continuous(limits = matdt[highlight=='high',range(eval(as.name(valueid)), na.rm = T)]) +
    # scale_alpha_discrete(range=c(1, 0.01)) +
    guides(colour=F,
           alpha=F)
  return(pl)
}

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

parseGI0 <- function(geneInfo0, groupcolumn='moduleColor', universe) {
  dt <- data.table(id=universe[!universe%in%geneInfo0[,id]],group=rep('background', times = length(universe[!universe%in%geneInfo0[,id]])))
  setnames(dt, 'group', groupcolumn)
  grouptable <- rbindlist(list(geneInfo0, dt))
  return(grouptable)
}

# https://rdrr.io/github/bryanhanson/ChemoSpecMarkeR/man/findElbow.html
# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve/2022348#2022348
findTOMCUT <- function(tom) {
  require(ChemoSpecMarkeR)
  m1 <- ecdf(tom)
  sqc <- seq(min(tom),max(tom),length.out = 1000)
  sqc <- sqc[2:(length(sqc)-1)]
  tomseq <- unlist(mclapply(sqc, m1, mc.cores = 10))
  tomseq1 <- tomseq[-length(tomseq)]
  stg <- lm(diff(tomseq)~tomseq1)
  upper <- seq(min(tom),max(tom),length.out = 1000)[findElbow(predict(stg, newdata = list(tomseq)))]
  sqc <- seq(min(tom),upper,length.out = 1000)
  sqc <- sqc[2:(length(sqc))]
  tomseq <- unlist(mclapply(sqc, m1, mc.cores = 10))
  # tomseq1 <- tomseq[-length(tomseq)]
  # stg <- lm(diff(tomseq)~tomseq1)
  # tomseq <- unlist(mclapply(seq(min(tom),upper,length.out = 1000), m1, mc.cores = 10))
  res <- seq(min(tom),upper,length.out = 1000)[findElbow(tomseq, plot = T)]
  return(res)
}

# require(data.table)
# geneInfo0 <- fread('res/20170831_wgcna/CSVs/wgcna.csv')
# universe <- readRDS('res/20170831_wgcna/RDS/universe.rds')
# grouptable <- parseGI0(geneInfo0, 'moduleColor', universe)
# func_ann <- readRDS('preproc/ath_ann.rds')

castfun <- function(identifier, keys, dt) {
  require(parallel)
  require(reshape2)
  res <- mclapply(dt[,unique(ONTOLOGY)],function(i) {
    setkeyv(dt, identifier)
    # TODO - needs to be amended using generateTermAssignment from jaccard.R
    tmp <- dt[.(keys),][!is.na(ID) & ONTOLOGY==i,list(eval(as.name(identifier)),ID,present)]
    setkey(tmp)
    tmp <- unique(tmp)
    tmp <- data.frame(acast(tmp, V1~ID, value.var = 'present', fill = F), check.names = F)
    tmp2 <- matrix(rep(F, times = ncol(tmp)*length(setdiff(keys, rownames(tmp)))), ncol = ncol(tmp))
    rownames(tmp2) <- setdiff(keys, rownames(tmp))
    colnames(tmp2) <- colnames(tmp)
    tmp <- rbind(tmp, tmp2)
    # apply(tmp, 2, function(i) factor(i, levels = unique(i), labels = unique(i), ordered = T))
    return(tmp)
  }, mc.cores = 10)
  names(res) <- dt[,unique(ONTOLOGY)]
  return(res)
}

quickGO <- function(grouptable,
                    gsubid = '_p[STY][0-9]*(_[0-9])?$',
                    subtypename='moduleColor',
                    func_ann = func_ann,
                    identifier = 'TAIR',
                    alternative="two.sided",
                    method.subtype = "none",
                    method.all = 'BH',
                    digits = 3,
                    cutoff = 0.05,
                    multimembersep = NULL,
                    subtypestoexclude = 'background',
                    go_rollup = func_ann$go_offspring,
                    ncores = 10) {
  grouptable[,id:=gsub(pattern = gsubid, replacement = '', id)]
  setkey(grouptable)
  grouptable <- unique(grouptable)
  if (nrow(grouptable)!=grouptable[,length(unique(id))]) {
    cat('Forced multimembersep to be "_" since identifiers are members of multiple groups!\n')
    multimembersep <- '_'
  }
  setkey(grouptable, id)
  grouptable <- grouptable[,list(tmpcol=paste0(eval(as.name(subtypename)), collapse = multimembersep)),by=id]
  setnames(grouptable, 'tmpcol', subtypename)
  tmp <- generateTermAssignment(ids = grouptable[,id], func_ann = func_ann, go_rollup = go_rollup, identifier = identifier)
  tmp <- unique(tmp[,list(eval(as.name(identifier)),PID)])
  setnames(tmp, 'V1', identifier)
  setkey(tmp, PID)
  setkey(func_ann$ontomap, ID)
  tmp[func_ann$ontomap,ONTOLOGY:=ONTOLOGY]
  tmp[,present:=T]
  tmp[is.na(ONTOLOGY),ONTOLOGY:='NA']
  setkey(tmp, ONTOLOGY)
  anntab2 <- pbmclapply(tmp[,unique(ONTOLOGY)], function(i) {
    tmp2 <- data.frame(acast(data = tmp[.(i),], formula = eval(as.name(identifier))~PID, value.var = 'present', fill = F), check.names = F)
    iddiff <- setdiff(grouptable[,id], rownames(tmp2))
    tmp2 <- rbind(tmp2, as.data.frame(matrix(F, nrow = length(iddiff), ncol = ncol(tmp2), dimnames = list(rownames=iddiff, colnames=colnames(tmp2)))))
  }, mc.cores = tmp[,uniqueN(ONTOLOGY)])
  names(anntab2) <- tmp[,unique(ONTOLOGY)]
  for (i in names(anntab2)) {
    anntab2[[i]][,subtypename] <- grouptable[.(rownames(anntab2[[i]])),eval(as.name(subtypename))] # HERE!
  }
  anntab <- anntab2
  # anntab <- castfun(identifier = identifier, keys = grouptable[,id], dt = func_ann$annotations)
  # for (i in names(anntab)) {
  #   anntab[[i]][,subtypename] <- grouptable[.(rownames(anntab[[i]])),eval(as.name(subtypename))] # HERE!
  # }
  res <- lapply(names(anntab), function(i) {
    cat(i, '\n')
    fisher.rbind(subtypename = subtypename,
                 categories = colnames(anntab[[i]])[!colnames(anntab[[i]])%in%subtypename],
                 func_ann = func_ann,
                 data = anntab[[i]],
                 alternative = alternative,
                 method.subtype = method.subtype,
                 method.all = method.all,
                 digits = digits,
                 cutoff = cutoff,
                 multimembersep = multimembersep,
                 subtypestoexclude = subtypestoexclude,
                 go_rollup = NULL,
                 ncores = ncores,
                 parseResult = T)
  })
  names(res) <- names(anntab)
  return(res)
}

parseGeneInfo0 <- function(geneInfo0, groupcolumn='moduleColor', universe) {
  setkeyv(geneInfo0, groupcolumn)
  tmp <- lapply(unique(geneInfo0[,eval(as.name(groupcolumn))]), function(i) as.factor(as.integer(universe%in%geneInfo0[.(i),id])))
  names(tmp) <- unique(geneInfo0[,eval(as.name(groupcolumn))])
  for (i in seq_along(tmp)) {
    names(tmp[[i]]) <- universe
  }
  return(tmp)
}

parseGOlist <- function(golist, universe) {
  require(parallel)
  require(topGO)
  require(annotate)
  geneid2go <- lapply(c('BP', 'MF', 'CC'), function(j) mclapply(universe, function(i) try(getOntology(golist[[i]], ontology = j), silent = T), mc.cores = 20))
  names(geneid2go) <- c('BP', 'MF', 'CC')
  for (i in seq_along(geneid2go)) {
    names(geneid2go[[i]]) <- universe
  }
  return(geneid2go)
}

parseGeneId2Go <- function(geneid2go) {
  require(parallel)
  allgo <- lapply(geneid2go, function(i) unique(unlist(i)))
  tmp <- lapply(seq_along(geneid2go), function(i) data.frame(do.call(rbind, mclapply(geneid2go[[i]], function(j) allgo[[i]]%in%j, mc.cores = 20))))
  for (i in seq_along(tmp)) {
    colnames(tmp[[i]]) <- allgo[[i]]
    tmp[[i]]$AGI <- rownames(tmp[[i]])
  }
  tmp <- Reduce(function(x, y) merge(x, y, by="AGI"), tmp)
  return(tmp)
}

doGO <- function(universe, geneInfo0, groupcolumn = 'moduleColor', nodesize, golist, cutoff = 0.05, outpath) {
  siggenelist <- parseGeneInfo0(geneInfo0, groupcolumn, universe)
  geneid2go <- parseGOlist(golist, universe)
  dir.create(file.path(outpath, 'RDS'), recursive = T, showWarnings = F)
  dir.create(file.path(outpath, 'CSVs', 'GO'), recursive = T, showWarnings = F)
  saveRDS(golist, file.path(outpath, 'RDS', 'golist.rds'), compress = T)
  gopath <- file.path(outpath, 'CSVs', 'GO')
  lapply(c('BP', 'MF', 'CC'), function(i) {
    lapply(names(siggenelist), function(j) {
      go <- try(new("topGOdata",
                    description = "First test",
                    ontology = i,
                    allGenes = siggenelist[[j]],
                    nodeSize = nodesize,
                    annot = annFUN.gene2GO,
                    gene2GO = geneid2go[[i]]), silent = T)
      if (class(go)!='try-error') {
        resultFisher <- runTest(go, algorithm = "classic", statistic = "fisher")
        allRes <- data.table(GenTable(go,
                                      classicFisher = resultFisher,
                                      orderBy = "classicFisher",
                                      ranksOf = "classicFisher",
                                      topNodes = length(resultFisher@score)))[classicFisher<=cutoff,]
        fwrite(allRes, file.path(gopath, paste0(i, '_', j, '.csv')), sep = ',', col.names = T, row.names = F) 
      }
    })
  })
}

doWGCNA <- function(mat1,
                    mat2=NULL,
                    setLabels=c('mat1'),
                    type = 'signed',
                    traits,
                    func_ann,
                    identifier = 'TAIR',
                    softPower = NULL,
                    MEDissThres = NULL,
                    agglomeration = 'average',
                    deepSplit = 2,
                    pamRespectsDendro = FALSE,
                    colord = NULL,
                    selectedTraits = colnames(traits),
                    plottedTraits = head(colnames(traits), 2),
                    outpath = getwd(),
                    scale = T,
                    scaleP = 0.95,
                    scalecut = 0,
                    sampleord = 'hclust',
                    excludeoutliers = T,
                    removenas = T,
                    minNOBS = 7,
                    filtermad = T,
                    quantcut = 0.75,
                    heatmaps = F,
                    writeentirenetwork = F,
                    automatedtomcut = F,
                    tomcut = 0.02,
                    minModuleSize = 30,
                    # nodesize = 10,
                    cutoff = 0.05,
                    gsubid = '_p[STY][0-9]*(_[0-9])?$') {
  cat(outpath, '\n')
  if (length(unique(setLabels))!=length(setLabels)) {
    stop('setLabels must be unique!')
  }
  
  if (!is.null(mat2) & !identical(rownames(mat1),rownames(mat2))) {
    stop('Rownames of mat1 and mat2 must be identical!')
  }
  dir.create(file.path(outpath, 'RDS'), recursive = T, showWarnings = F)
  # dir.create(file.path(outpath, 'RDATA'), recursive = T, showWarnings = F)
  dir.create(file.path(outpath, 'PLOTS'), recursive = T, showWarnings = F)
  dir.create(file.path(outpath, 'CSVs', 'GO'), recursive = T, showWarnings = F)
  dir.create(file.path(outpath, 'CSVs', 'Modules'), recursive = T, showWarnings = F)
  # source("/media/msdata5/users_files/Martin Frejno/phd/data/ms/24_CRC64_1and2and3_120ppm_1.4.1.2_peakproperties/analysis/fisher.R")
  source('/media/kusterlab/internal_projects/active/CRC65/data/24_CRC64_1and2and3_120ppm_1.4.1.2_peakproperties/analysis/colors.R')
  require(WGCNA)
  enableWGCNAThreads(nThreads = 16)
  require(data.table)
  require(pheatmap)
  require(arrayQualityMetrics)
  require(limma)
  # require(topGO)
  require(annotate)
  require(parallel)
  require(dendsort)
  options(stringsAsFactors = F)
  
  # Check last run
  
  lastmat1 <- NULL
  lastmat2 <- NULL
  lasttraits <- NULL
  # lastgolist <- NULL
  lastfunc_ann <- NULL
  # lastselectedTraits <- NULL
  # lastplottedTraits <- NULL
  lastcall <- NULL
  # lastsoftPower <- NULL
  # lastMEDissThres <- NULL
  # lastcolord <- NULL
  # lastscale <- NULL
  
  options(warn = -1)
  try(lastmat1 <- readRDS(file.path(outpath, 'RDS', 'mat.rds')), silent = T)
  try(lastmat2 <- readRDS(file.path(outpath, 'RDS', 'mat.rds')), silent = T)
  try(lasttraits <- readRDS(file.path(outpath, 'RDS', 'traits.rds')), silent = T)
  # try(lastgolist <- readRDS(file.path(outpath, 'RDS', 'golist.rds')), silent = T)
  try(lastfunc_ann <- readRDS(file.path(outpath, 'RDS', 'func_ann.rds')), silent = T)
  # try(lastselectedTraits <- readRDS(file.path(outpath, 'RDS', 'selectedTraits.rds')), silent = T)
  # try(lastplottedTraits <- readRDS(file.path(outpath, 'RDS', 'plottedTraits.rds')), silent = T)
  # try(lastsoftPower <- readRDS(file.path(outpath, 'RDS', 'softPower.rds')), silent = T)
  # try(lastMEDissThres <- readRDS(file.path(outpath, 'RDS', 'MEDissThres.rds')), silent = T)
  # try(lastcolord <- readRDS(file.path(outpath, 'RDS', 'colord.rds')), silent = T)
  # try(lastscale <- readRDS(file.path(outpath, 'RDS', 'scale.rds')), silent = T)
  try(lastcall <- readRDS(file.path(outpath, 'RDS', 'arguments.rds')), silent = T)
  options(warn = 0)
  
  # print(lastcall)
  
  # cat('HERE\n')
  # print(as.list(match.call()))
  
  funcall <- as.list(match.call())[5:length(as.list(match.call()))]
  # print(funcall)
  
  call <- unlist(lapply(1:length(funcall), function(i) !identical(eval(funcall[[i]]), eval(lastcall[[i]]))))
  names(call) <- names(funcall)
  
  runnew <- c(mat1=!identical(lastmat1, mat1),
              mat2=!identical(lastmat2, mat2),
              traits=!identical(lasttraits, traits),
              # golist=!identical(lastgolist, golist),
              func_ann=!identical(lastfunc_ann, func_ann),
              # softPower=!identical(lastsoftPower, softPower),
              # MEDissThres=!identical(lastMEDissThres,MEDissThres),
              # colord=!identical(lastcolord,colord),
              # selectedTraits=!identical(lastselectedTraits, selectedTraits),
              # plottedTraits=!identical(lastplottedTraits, plottedTraits),
              # scale=!identical(lastscale, scale),
              call)
  rm(lasttraits)
  rm(lastmat1)
  rm(lastmat2)
  gc()
  
  if(is.null(lastcall) || length(runnew)==0) {
    tmp <- rep(T, length(funcall))
    names(tmp) <- names(funcall)
    runnew <- c(mat1=T,
                mat2=T,
                traits=T,
                # golist=T,
                func_ann=T,
                # selectedTraits=T,
                # plottedTraits=T,
                # scale=T,
                # softPower=T,
                # MEDissThres=T,
                # colord=T,
                tmp)
  }
  print(runnew)
  
  saveRDS(mat1, file = file.path(outpath, 'RDS', 'mat1.rds'), compress = T)
  saveRDS(mat2, file = file.path(outpath, 'RDS', 'mat2.rds'), compress = T)
  saveRDS(traits, file = file.path(outpath, 'RDS', 'traits.rds'), compress = T)
  saveRDS(func_ann, file = file.path(outpath, 'RDS', 'func_ann.rds'), compress = T)
  # saveRDS(selectedTraits, file = file.path(outpath, 'RDS', 'selectedTraits.rds'), compress = T)
  # saveRDS(plottedTraits, file = file.path(outpath, 'RDS', 'plottedTraits.rds'), compress = T)
  # saveRDS(scale, file = file.path(outpath, 'RDS', 'scale.rds'), compress = T)
  saveRDS(funcall, file = file.path(outpath, 'RDS', 'arguments.rds'), compress = T)
  # saveRDS(colord, file.path(outpath, 'RDS', 'colord.rds'), compress = T)
  
  universe <- rownames(mat1)
  saveRDS(universe, file.path(outpath, 'RDS', 'universe.rds'), compress = T)
  
  
  # Form multi-set expression data: columns starting from 9 contain actual expression data.
  # We work with two sets:
  if (!is.null(mat2)) {
    nSets = 2
  } else {
    nSets = 1
  }
  multiExpr = vector(mode = "list", length = nSets)
  # For easier labeling of plots, create a vector holding descriptive names of the two sets.
  multiExpr[[1]] = list(data = as.data.frame(t(mat1)))
  if (!is.null(mat2)) {
    multiExpr[[2]] = list(data = as.data.frame(t(mat2)))
  }
  shortLabels = setLabels
  names(multiExpr) <- setLabels
  
  exprSize = checkSets(multiExpr)
  
  
  # Prefilter matrix based on mad
  if (any(runnew['mat1'],runnew['mat2'],runnew['quantcut'],runnew['filtermad']) & filtermad) {
    cat('Reduce matrix!\n')
    mads <- lapply(multiExpr, function(j) apply(j$data, 2, function(i) mad(i, na.rm = T)))
    qcut <- lapply(mads, function(i) names(which(i>=quantile(x = i, probs = quantcut, na.rm = T))))
    gtk <- Reduce(intersect, qcut)
    for (i in 1:nSets) {
      multiExpr[[i]]$data <- multiExpr[[i]]$data[,gtk]
    }
  }
  
  exprSize = checkSets(multiExpr)
  
  if (runnew['removenas'] & removenas) {
    gsg <- goodSamplesGenesMS(multiExpr, verbose = 3)
    if (!gsg$allOK) {
      # Optionally, print the gene and sample names that were removed:
      # if (sum(!gsg$goodGenes)>0) 
      #   printFlush(paste("Removing genes:", paste(names(mat)[!gsg$goodGenes], collapse = ", ")));
      # if (sum(!gsg$goodSamples)>0) 
      #   printFlush(paste("Removing samples:", paste(rownames(mat)[!gsg$goodSamples], collapse = ", ")));
      # Remove the offending genes and samples from the data:
      for (i in 1:exprSize$nSets) {
        multiExpr[[i]]$data <- multiExpr[[i]]$data[gsg$goodSamples[[i]],gsg$goodGenes]
      }
    }
    saveRDS(gsg, file.path(outpath, 'RDS', 'gsg.rds'), compress = T) 
  }
  
  # remove outliers
  if (any(runnew['mat1'],runnew['mat2'],runnew['quantcut'],runnew['excludeoutliers']) & excludeoutliers) {
    dir.create(file.path(outpath, 'AQC'), recursive = T, showWarnings = F)
    for (i in 1:exprSize$nSets) {
      exp_mat <- new(Class = 'ExpressionSet', exp = as.matrix(t(multiExpr[[i]]$data)))
      options(warn = -1)
      qc <- arrayQualityMetrics(expressionset = exp_mat, outdir = file.path(outpath, 'AQC', setLabels[i]), force = T, do.logtransform = F)
      try(dev.off(), silent = T)
      options(warn = 0)
      saveRDS(qc, file.path(outpath, 'RDS', paste0(setLabels[i],'_qc.rds')), compress = T)
      out <- rowSums(qc$arrayTable[,3:5]=="x")!=0
      multiExpr[[i]]$data <- multiExpr[[i]]$data[!out,]
    }
    saveRDS(multiExpr, file.path(outpath, 'RDS', 'multiExpr_nooutliers.rds'), compress = T)
  } else {
    saveRDS(multiExpr, file.path(outpath, 'RDS', 'multiExpr_nooutliers.rds'), compress = T)
  }
  
  exprSize = checkSets(multiExpr)
  
  multiExpr <- readRDS(file.path(outpath, 'RDS', 'multiExpr_nooutliers.rds'))
  universe <- readRDS(file.path(outpath, 'RDS', 'universe.rds'))
  
  if (runnew['scale']) {
    if (scale) {
      for (i in 1:nSets) {
        multiExpr[[i]]$data <- scale(multiExpr[[i]]$data)
      }
    } else {
      multiExpr <- multiExpr
    }
    saveRDS(multiExpr, file.path(outpath, 'RDS', 'multiExpr_afterscale.rds'), compress = T)
  }
  
  candp <- NULL
  try(candp <- readRDS(file.path(outpath, 'RDS', 'corandp.rds')), silent = T)
  
  if (is.null(candp)) {
    cat('Calculate all pairwise correlations!\n')
    candp <- lapply(multiExpr, function(i) {
        corandp(i$data, use = "pairwise.complete.obs", alternative = "two.sided", nThreads = 16)
      })
    cat('Write correlation matrix to disk!\n')
    saveRDS(candp, file.path(outpath, 'RDS', 'corandp.rds'), compress = T)
  }
  
  if (minNOBS>0) {
    cat('Mask pairs with low overlap!\n')
    mask <- Reduce(`|`, lapply(candp, function(i) i$nObs<minNOBS))
    for (i in names(candp)) {
      candp[[i]]$cor[mask] <- NA
    }
  }
  
  simfun <- function(cor, type) {
    type <- match.arg(type, c('unsigned', 'signed'))
    switch(type,
           'unsigned' = abs(cor),
           'signed' = (0.5*(1+cor)))
  }
  
  simlist <- lapply(candp, function(i) simfun(cor = i$cor, type = type))
  # names(simlist) <- setLabels
  
  rm(candp)
  gc()
  
  if (is.null(softPower)) {
    # source('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/src/pickSoftThresholdFIX.R')
    # Choose a set of soft-thresholding powers
    powers = c(c(1:10), seq(from = 12, to=30, by=2))
    # Call the network topology analysis function
    powerTables <- list()
    cat('Calculate power tables!\n')
    for (set in 1:nSets) {
      powerTables[[set]] = list(data = pickSoftThreshold(simlist[[set]], dataIsExpr = F, powerVector=powers, verbose = 2, networkType = type, blockSize = min(ceiling(c(2000, nrow(simlist[[set]])/5))))[[2]])
    }
    
    # Plot the results:
    colors = c("black", "red")
    # Will plot these columns of the returned scale free analysis tables
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
                 "Max connectivity");
    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4);
    for (set in 1:nSets)
    {
      for (col in 1:length(plotCols))
      {
        ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
      }
    }
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    sizeGrWindow(8, 6)
    png(file.path(outpath, 'PLOTS', 'scale_free_topology.png'), width = (210-32.5*2)*2, height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
    # pdf(file = "Plots/scaleFreeAnalysis.pdf", wi = 8, he = 6)
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;
    for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
      if (set==1)
      {
        plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
             main = colNames[col]);
        addGrid();
      }
      if (col==1)
      {
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
             labels=powers,cex=cex1,col=colors[set]);
      } else
        text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
             labels=powers,cex=cex1,col=colors[set]);
      if (col==1)
      {
        legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
      } else
        legend("topright", legend = setLabels, col = colors, pch = 20) ;
    }
    dev.off();
    
  }
  
  multiExpr <- readRDS(file.path(outpath, 'RDS', 'multiExpr_afterscale.rds'))
  
  if (runnew['traits']) {
    # traits <- traits[rownames(mat),]
    # Re-cluster samples
    sampleCors = list()
    for (set in 1:exprSize$nSets)
    {
      sampleCors[[set]] = cor(t(multiExpr[[set]]$data), use = 'p')
      # sampleCors[[set]] = hclust(as.dist(1-sampleCors[[set]]), method = agglomeration)
      # sampleCors[[set]] = dendsort(sampleCors[[set]], type = 'min')
    }
    sampleTrees = list()
    for (set in seq_along(sampleCors))
    {
      # sampleTrees[[set]] = cor(t(multiExpr[[set]]$data), use = 'p')
      sampleTrees[[set]] = hclust(as.dist(1-sampleCors[[set]]), method = agglomeration)
      sampleTrees[[set]] = dendsort(sampleTrees[[set]], type = 'min')
    }
    names(sampleTrees) <- setLabels
    if (exprSize$nSets>1) {
      sampleTrees[['consTOM']] <- dendsort(hclust(as.dist(1-cor(sampleCors[[1]], sampleCors[[2]])), method = agglomeration), type = 'min')
    }
    # sampleTree2 = hclust(as.dist(1-cor(t(mat), use = 'p')), method = agglomeration)
    # sampleTree2 <- dendsort(sampleTree2, type = 'min')
    # Convert traits to a color representation: white means low, red means high, grey means missing entry
    traitColors = labels2colors(traits)
    # Plot the sample dendrogram and the colors underneath.
    png(file.path(outpath, 'PLOTS', 'sampleclustering.png'), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
    par(mfrow=c(2,1))
    par(mar = c(0, 4, 2, 0))
    for (set in 1:nSets)
      plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
           xlab="", sub="", cex = 0.7);
    dev.off()
    
    saveRDS(sampleTrees, file.path(outpath, 'RDS', 'sampleTrees.rds'), compress = T)
  }
  
  if (any(runnew['softPower'],runnew['minModuleSize'])) {
    options(warn = -1)
    while (is.null(softPower) || is.na(softPower) || softPower%%1!=0 || softPower>30) {
      softPower <- as.numeric(readline(prompt = 'Pick soft threshold power:\n'))
      if (!is.na(softPower) & softPower%%1==0 & softPower<=30) {
      } else {
        cat('Please select an integer between 1 and 30!')
      }
    }
    options(warn = 0)
    cat('Build network!\n')
    
    saveRDS(softPower, file.path(outpath, 'RDS', 'softPower.rds'), compress = T)
    
    exprSize <- checkSets(multiExpr)
    
    cat('Calculate adjacency!\n')
    # adjlist <- lapply(multiExpr, function(i) adjacency(i$data, power = softPower, distOptions = list(method = agglomeration)))

    adjfun <- function(sim, softPower) {
      sim^softPower
    }
    adjlist <- lapply(simlist, function(i) adjfun(sim = i, softPower = softPower))
    names(adjlist) <- setLabels
    
    # if (minNOBS>0) {
    #   cat('Mask pairs with low overlap!\n')
    #   mask <- Reduce(`|`, lapply(candp, function(i) i$nObs<minNOBS))
    #   for (i in names(adjlist)) {
    #     adjlist[[i]][mask] <- NA
    #   }
    # }
    
    
    cat('Calculate TOM!\n')
    TOMlist <- lapply(adjlist, TOMsimilarity)
    names(TOMlist) <- setLabels
    
    
    # adjacency = adjacency(mat, power = softPower, distOptions = list(method = agglomeration))
    # Initialize an appropriate array to hold the adjacencies
    # adjacencies = array(0, dim = c(exprSize$nSets, exprSize$nGenes, exprSize$nGenes));
    # # Calculate adjacencies in each individual data set
    # for (set in 1:nSets) adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower
    # dimnames(adjacencies) <- list(setLabels, colnames(multiExpr[[1]]$data),colnames(multiExpr[[1]]$data))
    # 
    # # Turn adjacency into topological overlap
    # # TOM = TOMsimilarity(adjacency)
    # # Initialize an appropriate array to hold the TOMs
    # TOM = array(0, dim = c(exprSize$nSets, exprSize$nGenes, exprSize$nGenes))
    # # Calculate TOMs in each individual data set
    # for (set in 1:nSets) TOM[set, , ] = TOMsimilarity(adjacencies[set, , ])
    # 
    # dimnames(TOM) <- list(setLabels, colnames(multiExpr[[1]]$data),colnames(multiExpr[[1]]$data))
    
    # saveRDS(adjlist, file.path(outpath, 'RDS', 'adjacencies.rds'), compress = T)
    saveRDS(TOMlist, file.path(outpath, 'RDS', 'TOMlist.rds'), compress = T)
    
    if (!is.null(mat2) & runnew['scaleP']) {
      cat('Scale TOM!\n')
      # Define the reference percentile
      scaleP = scaleP
      # Set RNG seed for reproducibility of sampling
      set.seed(12345)
      # Sample sufficiently large number of TOM entries
      nSamples = as.integer(1/(1-scaleP) * 1000);
      # Choose the sampled TOM entries
      scaleSample = sample(exprSize$nGenes*(exprSize$nGenes-1)/2, size = nSamples)
      TOMScalingSamples = list();
      # These are TOM values at reference percentile
      scaleQuant = rep(1, exprSize$nSets)
      # Scaling powers to equalize reference TOM values
      scalePowers = rep(1, exprSize$nSets)
      # Loop over sets
      for (set in 1:exprSize$nSets)
      {
        # Select the sampled TOM entries
        TOMScalingSamples[[set]] = as.dist(TOMlist[[set]])[scaleSample]
        # Calculate the 95th percentile
        scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                                   probs = scaleP, type = 8);
        # Scale the male TOM
        if (set>1)
        {
          scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
          TOMlist[[set]] = TOMlist[[set]]^scalePowers[set];
        }
      }
      
      
      #=====================================================================================
      #
      #  Code chunk 6
      #
      #=====================================================================================
      
      
      # For plotting, also scale the sampled TOM entries
      scaledTOMSamples = list();
      for (set in 1:nSets)
        scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
      # Open a suitably sized graphics window
      sizeGrWindow(6,6)
      png(file.path(outpath, 'PLOTS', 'TOMScaling-QQPlot.png'), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
      # pdf(file = "Plots/TOMScaling-QQPlot.pdf", wi = 6, he = 6);
      # qq plot of the unscaled samples
      qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                          xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                          main = "Q-Q plot of TOM", pch = 20)
      # qq plot of the scaled samples
      qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
      points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
      abline(a=0, b=1, col = "blue")
      legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
      dev.off();
      TOMlist[['consTOM']] = pmin(TOMlist[[1]], TOMlist[[2]])
      saveRDS(TOMlist, file.path(outpath, 'RDS', 'TOMlist.rds'), compress = T)
    }
    
    if (!exists('TOMlist')) {
      TOMlist <- readRDS(file.path(outpath, 'RDS', 'TOMlist.rds'))
    }
    
    # bla <- array(sapply(1:dim(TOM)[1], function(i) 1-TOM[i,,]), dim = dim(TOM))
    # dissTOM = 1-TOM
    
    dissTOMlist <- lapply(TOMlist, function(i) 1-i)
    names(dissTOMlist) <- names(TOMlist)
    for (i in seq_along(dissTOMlist)) {
      rownames(dissTOMlist[[i]]) <- colnames(multiExpr[[1]]$data)
      colnames(dissTOMlist[[i]]) <- colnames(multiExpr[[1]]$data)
    }
    
    # Call the hierarchical clustering function
    cat('Calculate geneTree!\n')
    geneTreelist = lapply(dissTOMlist, function(i) hclust(as.dist(i), method = agglomeration))
    # cat('Sort geneTree!\n')
    # geneTreelist = lapply(geneTreelist, function(i) dendsort(i, type = 'min'))
    
    # Plot the resulting clustering tree (dendrogram)
    # plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
    #      labels = FALSE, hang = 0.04)
    
    # We like large modules, so we set the minimum module size relatively high:
    # minModuleSize = 30
    
    # Module identification using dynamic tree cut:
    cat('Identify modules!\n')
    dynamicModslist = lapply(seq_along(geneTreelist), function(i) cutreeDynamic(dendro = geneTreelist[[i]], distM = dissTOMlist[[i]],
                                                                                deepSplit = deepSplit, pamRespectsDendro = pamRespectsDendro,
                                                                                minClusterSize = minModuleSize))
    names(dynamicModslist) <- names(geneTreelist)
    
    # table(dynamicMods)
    
    # Convert numeric lables into colors
    dynamicColorslist = lapply(dynamicModslist, labels2colors)
    # table(dynamicColors)
    # Plot the dendrogram and colors underneath
    lapply(names(dynamicColorslist), function(i) {
      png(file.path(outpath, 'PLOTS', paste0('geneTree_pre_merge_', i, '.png')), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
      plotDendroAndColors(geneTreelist[[i]], dynamicColorslist[[i]], "Dynamic Tree Cut",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = paste0("Gene dendrogram and\nmodule colors - ", i))
      dev.off()
    })
    
    if (!is.null(mat2)) {
      # Calculate module eigengenes
      unmergedMEs = multiSetMEsFIX(multiExpr, colors = NULL, universalColors = dynamicColorslist$consTOM, impute = F, softPower = softPower, scale = !scale)
      names(unmergedMEs) <- setLabels
      # Calculate consensus dissimilarity of consensus module eigengenes
      consMEDiss = consensusMEDissimilarity(unmergedMEs)
      # Cluster consensus modules
      METree = hclust(as.dist(consMEDiss), method = agglomeration)
      METree = dendsort(METree, type = 'min')
      # Plot the result
      # sizeGrWindow(7,6)
      # par(mfrow = c(1,1))
      png(file.path(outpath, 'PLOTS', 'ME_pre_merge_cons.png'), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
      plot(METree, main = "Consensus clustering of consensus module eigengenes",
           xlab = "", sub = "")
      dev.off()
      # abline(h=0.25, col = "red")
      saveRDS(METree, file.path(outpath, 'RDS', 'METree.rds'), compress = T)
    }
    # Calculate eigengenes
    MEList = lapply(setLabels, function(i) moduleEigengenes(multiExpr[[i]]$data, colors = dynamicColorslist[[i]], softPower = softPower, impute = F, scale = !scale))
    names(MEList) <- setLabels
    MEs = lapply(MEList, function(i) i$eigengenes)
    # Calculate dissimilarity of module eigengenes
    MEDisslist = lapply(MEs, function(i) 1-cor(i, use = 'p'))
    # Cluster module eigengenes
    METreelist = lapply(MEDisslist, function(i) hclust(as.dist(i), method = agglomeration))
    # METree %>% as.dendrogram() %>% ladderize() %>% plot()
    METreelist <- lapply(METreelist, dendsort, type = 'min')
    # Plot the result
    lapply(names(METreelist), function(i) {
      png(file.path(outpath, 'PLOTS', paste0('ME_pre_merge_', i, '.png')), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
      plot(METreelist[[i]], main = paste0("Clustering of module eigengenes\n", i),
           xlab = "", sub = "")
      dev.off()
    })
  
    saveRDS(METreelist, file.path(outpath, 'RDS', 'METreelist.rds'), compress = T)
    saveRDS(dynamicColorslist, file.path(outpath, 'RDS', 'dynamicColors.rds'), compress = T)
    saveRDS(geneTreelist, file.path(outpath, 'RDS', 'geneTree.rds'), compress = T)
    # # rm(geneTreelist)
    # # rm(TOMlist)
    # # rm(dissTOMlist)
    # # rm(adjlist)
    # gc()
  }
  
  
  
  if (runnew['MEDissThres']) {
    options(warn = -1)
    while (is.null(MEDissThres) || is.na(MEDissThres) || (MEDissThres>1 | MEDissThres<0)) {
      MEDissThres <- as.numeric(readline(prompt = 'Pick a cutoff to merge modules between 0 and 1:\n'))
      if (!is.na(MEDissThres) & (MEDissThres<=1 & MEDissThres>=0)) {
      } else {
        cat('Please select a cutoff between 0 and 1!')
      }
    }
    options(warn = 0)
    cat('Merge modules!\n')
    
    saveRDS(MEDissThres, file.path(outpath, 'RDS', 'MEDissThres.rds'), compress = T)
    
    multiExpr <- readRDS(file.path(outpath, 'RDS', 'multiExpr_afterscale.rds'))
    METreelist <- readRDS(file.path(outpath, 'RDS', 'METreelist.rds'))
    dynamicColorslist <- readRDS(file.path(outpath, 'RDS', 'dynamicColors.rds'))
    geneTreelist <- readRDS(file.path(outpath, 'RDS', 'geneTree.rds'))
    
    if (!is.null(mat2)) {
      METree <- readRDS(file.path(outpath, 'RDS', 'METree.rds'))
      png(file.path(outpath, 'PLOTS', 'ME_post_merge_cons.png'), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
      plot(METree, main = "Consensus clustering of consensus module eigengenes",
           xlab = "", sub = "")
      abline(h=MEDissThres, col = "red")
      dev.off()
      # abline(h=0.25, col = "red")
    }
    # Plot the result
    lapply(names(METreelist), function(i) {
      png(file.path(outpath, 'PLOTS', paste0('ME_post_merge_', i, '.png')), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
      plot(METreelist[[i]], main = paste0("Clustering of module eigengenes\n", i),
           xlab = "", sub = "")
      abline(h=MEDissThres, col = "red")
      dev.off()
    })
    
    genmatlist <- function(multiExpr) {
      # reslist <- vector(mode = 'list', length = length(multiExpr)+1)
      reslist <- list()
      for (i in names(multiExpr)) {
        reslist[[i]] <- multiExpr[[i]]$data
      }
      if (length(multiExpr)>1) {
        reslist[['consTOM']] <- multiExpr
      }
      return(reslist)
    }
    
    matlist <- genmatlist(multiExpr)
    names(matlist) <- names(dynamicColorslist)
    
    # Call an automatic merging function
    mergelist = lapply(names(dynamicColorslist), function(i) mergeCloseModules(matlist[[i]], dynamicColorslist[[i]], cutHeight = MEDissThres, verbose = 3, impute = F))
    names(mergelist) <- names(dynamicColorslist)
    # merge = mergeCloseModules(mat, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    # The merged module colors
    mergedColors = lapply(mergelist, function(i) i$colors)
    # Eigengenes of the new merged modules:
    mergedMEslist = lapply(mergelist, function(i)  i$newMEs)
    
    # MEs = lapply(MEList, function(i) i$eigengenes)
    # # Calculate dissimilarity of module eigengenes
    # MEDisslist = lapply(MEs, function(i) 1-cor(i, use = 'p'))
    # # Cluster module eigengenes
    # METreelist = lapply(MEDisslist, function(i) hclust(as.dist(i), method = agglomeration))
    # # METree %>% as.dendrogram() %>% ladderize() %>% plot()
    # METreelist <- lapply(METreelist, dendsort, type = 'min')
    # # Plot the result
    # lapply(names(METreelist), function(i) {
    #   png(file.path(outpath, 'PLOTS', paste0('ME_pre_merge_', i, '.png')), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
    #   plot(METreelist[[i]], main = paste0("Clustering of module eigengenes\n", i),
    #        xlab = "", sub = "")
    #   dev.off()
    # })
    
    # Calculate dissimilarity of module eigengenes
    MEDisslist = lapply(mergedMEslist, function(i) if(class(i)=='list') lapply(i, function(j) 1-cor(j$data, use = 'p')) else 1-cor(i, use = 'p'))
    # MEDisslist = lapply(mergedMEslist, function(i) 1-cor(i, use = 'p'))
    # Cluster module eigengenes
    METreelist = lapply(MEDisslist, function(i) if(class(i)=='list') lapply(i, function(j) hclust(as.dist(j), method = agglomeration)) else hclust(as.dist(i), method = agglomeration))
    # METree %>% as.dendrogram() %>% ladderize() %>% plot()
    METreelist <- lapply(METreelist, function(i) if(class(i)=='list') lapply(i, function(j) dendsort(j, type = 'min')) else (dendsort(i, type = 'min')))
    # Plot the result
    lapply(names(METreelist), function(i) {
      if (class(METreelist[[i]]) == 'list') {
        lapply(names(METreelist[[i]]), function(j) {
          pdf(file.path(outpath, 'PLOTS', paste0('ME_post_merge_new_modules_', i, '_', j, '.pdf')),pointsize = 10)
          plot(METreelist[[i]][[j]], main = paste0("Clustering of module eigengenes\n", i, '_', j),
               xlab = "", sub = "")
          dev.off() 
        })
      } else {
        pdf(file.path(outpath, 'PLOTS', paste0('ME_post_merge_new_modules_', i, '.pdf')),pointsize = 10)
        plot(METreelist[[i]], main = paste0("Clustering of module eigengenes\n", i),
             xlab = "", sub = "")
        dev.off()  
      }
    })
    require(reshape2)
    MEDissCSV <- rbindlist(lapply(names(MEDisslist), function(i) {
      if (class(MEDisslist[[i]])=='list') {
        rbindlist(lapply(names(MEDisslist[[i]]), function(j) {
          tmp <- data.table(melt(MEDisslist[[i]][[j]]))
          setnames(tmp, c('ME1', 'ME2', 'Dissimilarity'))
          tmp[,Dataset:=paste0(i, '_', j)]
          return(tmp)
        }))
      } else {
        tmp <- data.table(melt(MEDisslist[[i]]))
        setnames(tmp, c('ME1', 'ME2', 'Dissimilarity'))
        tmp[,Dataset:=i]
        return(tmp) 
      }
    }))
    
    fwrite(MEDissCSV, file.path(outpath, 'CSVs', 'CombinedMEDissimilarity.csv'), sep = ',', col.names = T, row.names = F)
    cols <- lapply(MEDisslist, function(i) if(class(i)=='list') lapply(i, function(j) gsub('ME', '', colnames(j))) else gsub('ME', '', colnames(i)))
    
    pdf(file = file.path(outpath, 'PLOTS', 'MergedModulesColorLegend.pdf'), onefile = T, pointsize = 8)
    par(mar=c(10.1,4.1,4.1,2.1))
    lapply(seq_along(cols), function(i) if(class(cols[[i]])=='list') lapply(seq_along(cols[[i]]), function(j) barplot(rep(1, length(cols[[i]][[j]])), col = cols[[i]][[j]], names.arg = cols[[i]][[j]], las = 2, yaxt = 'none', main = paste0('Legend for ', names(cols)[i], '_', names(cols[[i]])[j]))) else barplot(rep(1, length(cols[[i]])), col = cols[[i]], names.arg = cols[[i]], las = 2, yaxt = 'none', main = paste0('Legend for ', names(cols)[i])))
    dev.off()
    
    # sizeGrWindow(12, 9)
    lapply(names(dynamicColorslist), function(i) {
      png(file = file.path(outpath, 'PLOTS', paste0("geneTree_post_merge_", i, '.png')), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
      plotDendroAndColors(geneTreelist[[i]], cbind(dynamicColorslist[[i]], mergedColors[[i]]),
                          c("Dynamic Tree Cut", "Merged dynamic"),
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = paste0("Gene dendrogram and\nmodule colors - ", i))
      dev.off() 
    })
    
    if (!is.null(mat2)) {
      consMEDiss = consensusMEDissimilarity(mergedMEslist$consTOM)
      # Cluster consensus modules
      METree = hclust(as.dist(consMEDiss), method = agglomeration)
      METree = dendsort(METree, type = 'min')
      saveRDS(METree, file.path(outpath, 'RDS', 'METree.rds'), compress = T)
    }
    
    saveRDS(mergedColors, file.path(outpath, 'RDS', 'mergedColors.rds'), compress = T)
    saveRDS(mergedMEslist, file.path(outpath, 'RDS', 'mergedMEs.rds'), compress = T)
    saveRDS(matlist, file.path(outpath, 'RDS', 'matlist.rds'), compress = T)
  }
  
  matlist <- readRDS(file.path(outpath, 'RDS', 'matlist.rds'))
  mergedColors <- readRDS(file.path(outpath, 'RDS', 'mergedColors.rds'))
  geneTreelist <- readRDS(file.path(outpath, 'RDS', 'geneTree.rds'))
  mergedMEslist <- readRDS(file.path(outpath, 'RDS', 'mergedMEs.rds'))
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(max(sapply(moduleColors, function(i) length(unique(i))))))
  moduleLabels = lapply(moduleColors, function(i) match(i, colorOrder)-1)
  MEslist = mergedMEslist
  
  # Save module colors and labels for use in subsequent parts
  # save(MEs, moduleLabels, moduleColors, geneTree, file = "after_network_construction.RData")
  
  # Define numbers of genes and samples
  multiExpr <- readRDS(file.path(outpath, 'RDS', 'multiExpr_afterscale.rds'))
  nGenes = as.list(rep(checkSets(multiExpr)$nGenes, times = length(matlist)))
  names(nGenes) <- names(matlist)
  nSamples = lapply(matlist, function(i) if(is.integer(nrow(i))) nrow(i) else checkSets(i)$nSamples)
  if (!is.null(mat2)) {
    names(nSamples$consTOM) <- setLabels
  }
  # Recalculate MEs with color labels
  # Calculate eigengenes
  MEs0list = lapply(setLabels, function(i) moduleEigengenes(matlist[[i]], colors = moduleColors[[i]], softPower = softPower, impute = F, scale = !scale)$eigengenes)
  names(MEs0list) <- setLabels
  if (!is.null(mat2)) {
    # Calculate module eigengenes
    mergedMEscons = lapply(multiSetMEsFIX(matlist$consTOM, colors = NULL, universalColors = mergedColors$consTOM, impute = F, softPower = softPower, scale = !scale), function(i) i$data)
    # Calculate consensus dissimilarity of consensus module eigengenes
    # mergedconsMEDiss = consensusMEDissimilarity(mergedMEscons)
    # Cluster consensus modules
    # METree = hclust(as.dist(consMEDiss), method = agglomeration);
    # METree = dendsort(METree, type = 'min')
    # # Plot the result
    # # sizeGrWindow(7,6)
    # # par(mfrow = c(1,1))
    # png(file.path(outpath, 'PLOTS', 'ME_pre_merge_cons.png'), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
    # plot(METree, main = "Consensus clustering of consensus module eigengenes",
    #      xlab = "", sub = "")
    # dev.off()
    # abline(h=0.25, col = "red")
    MEs0list[['consTOM']] <- mergedMEscons
  }
  
  if (!is.null(mat2)) {
    MET <- vector(mode = 'list', length = length(matlist))
    names(MET) <- names(matlist)
    for (i in names(matlist)) {
      if (class(matlist[[i]])=='list') {
          MET[[i]] <- multiSetMEsFIX(matlist[[i]], colors = NULL, universalColors = mergedColors[[i]], impute = F, softPower = softPower, scale = !scale)
      } else {
        MET[[i]] <- moduleEigengenes(matlist[[i]], colors = moduleColors[[i]], softPower = softPower, impute = F, scale = !scale)
      }
    }
    # ordMET <- consensusOrderMEs(MET$consTOM)
    METree <- readRDS(file.path(outpath, 'RDS', 'METree.rds'))
    ordMET <- orderMEs(MET$consTOM, order = METree$order)
    
    png(file.path(outpath, 'PLOTS', paste0('EigengeneNetworks.png')), width = (210-32.5*2)*1.5, height = (210-32.5*2)*2, units = "mm", res = 300, pointsize = 6)
    plotEigengeneNetworks(multiME = ordMET,
                          excludeGrey = F,
                          plotAdjacency = F,
                          setLabels = names(ordMET),
                          marDendro = c(0,2,2,1),
                          marHeatmap = c(10,10,1,1),
                          zlimPreservation = c(0.5, 1),
                          xLabelsAngle = 90,
                          barplotErrors = T,
                          colorLabels = F,
                          coloredBarplot = T,
                          xLabels = colnames(ordMET$RNA$data),
                          yLabels = colnames(ordMET$RNA$data))
    dev.off()
  }
  
  # MEs0list = lapply(names(matlist), function(i) moduleEigengenes(matlist[[i]], moduleColors[[i]])$eigengenes)
  MEslist = lapply(MEs0list, function(i) if(class(i)=='list') lapply(i, orderMEs) else orderMEs(i))
  saveRDS(MEslist, file.path(outpath, 'RDS', 'mergedMEs.rds'))
  MEorderlist <- lapply(MEs0list, function(i) if(class(i)=='list') lapply(i, colnames) else colnames(i))
  saveRDS(MEorderlist, file.path(outpath, 'RDS', 'MEorder.rds'), compress = T)
  MEslistnonas = lapply(MEslist, function(i) if(class(i)=='list') lapply(i, function(j) j[,which(!is.nan(colSums(j)))]) else i[,which(!is.nan(colSums(i)))])
  iftest <- identical(lapply(MEslist, function(i) if(class(i)=='list') lapply(i, ncol) else ncol(i)),lapply(MEslistnonas, function(i) if(class(i)=='list') lapply(i, ncol) else ncol(i)))
  if (!iftest) {
    cat('Removed some MEs due to invalid values!')
  }
  
  if (!is.null(mat2)) {
    cat('Relating consensus modules to set-specific modules!\n')
    moduleColors <- readRDS(file.path(outpath, 'RDS', 'mergedColors.rds'))
    MEslist <- readRDS(file.path(outpath, 'RDS', 'mergedMEs.rds'))
    modlabs <- lapply(MEslist, function(i) if(class(i)=='list') lapply(i, function(j) substring(names(j), 3)) else substring(names(i), 3))
    countmod <- lapply(modlabs, function(i) if(class(i)=='list') lapply(i, length) else length(i))
    require(parallel)
    ptablelist <- mclapply(setLabels, function(set)
    {
      cat(set, '\n')
      reslist <- vector(mode = 'list', length = 4L)
      names(reslist) <- c('pTable', 'CountTbl', 'setModTotals', 'consModTotals')
      reslist[['pTable']] = matrix(0, nrow = countmod[[set]], ncol = countmod[['consTOM']][[set]])
      rownames(reslist[['pTable']]) <- modlabs[[set]]
      colnames(reslist[['pTable']]) <- modlabs[['consTOM']][[set]]
      reslist[['CountTbl']] = matrix(0, nrow = countmod[[set]], ncol = countmod[['consTOM']][[set]])
      rownames(reslist[['CountTbl']]) <- modlabs[[set]]
      colnames(reslist[['CountTbl']]) <- modlabs[['consTOM']][[set]]
      for (setmod in 1:countmod[[set]])
      {
        for (cmod in 1:countmod[['consTOM']][[set]])
        {
          setMembers = (moduleColors[[set]] == modlabs[[set]][[setmod]])
          consMembers = (moduleColors[["consTOM"]] == modlabs[["consTOM"]][[set]][[cmod]]);
          reslist[['pTable']][setmod, cmod] = -log10(fisher.test(setMembers, consMembers, alternative = "greater")$p.value)
          # reslist[['pTable']][is.infinite(reslist[['pTable']])] = 1.3*max(reslist[['pTable']][is.finite(reslist[['pTable']])])
          # reslist[['pTable']][reslist[['pTable']]>50 ] = 50
          reslist[['CountTbl']][setmod, cmod] = sum(setMembers & consMembers)
        }
      }
      reslist[['setModTotals']] <- apply(reslist[['CountTbl']], 1, sum)
      reslist[['consModTotals']] <- apply(reslist[['CountTbl']], 2, sum)
      return(reslist)
    }, mc.cores = length(setLabels)
    )
    names(ptablelist) <- setLabels
    require(reshape2)
    require(data.table)
    ptablist <- lapply(ptablelist, function(i) data.table(Reduce(merge, list(melt(i$pTable, value.name = 'p value'), melt(i$CountTbl, value.name = 'Overlap'))))[,`p value`:=10^(-`p value`)])
    lapply(names(ptablist), function(i) setnames(ptablist[[i]], c('Var1', 'Var2'), c(i, 'consTOM')))
    lapply(names(ptablist), function(i) fwrite(ptablist[[i]], file = file.path(outpath, 'CSVs', paste0('Module_Correspondence_', i, '_consTOM.csv')), col.names = T, row.names = F, sep = ','))
    for (i in names(ptablelist)) {
      ptablelist[[i]][['pTable']][is.infinite(ptablelist[[i]][['pTable']])] = 1.3*max(ptablelist[[i]][['pTable']][is.finite(ptablelist[[i]][['pTable']])])
      ptablelist[[i]][['pTable']][ptablelist[[i]][['pTable']]>50 ] = 50
    }
    lapply(names(ptablelist), function(i) {
      png(file.path(outpath, 'PLOTS', paste0('Correspondence_heatmap_', i, '.png')), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
      par(mfrow=c(1,1))
      par(cex = 0.5)
      par(mar=c(12, 12, 2.7, 1)+0.3)
      labeledHeatmap(Matrix = ptablelist[[i]]$pTable[order(paste(i, " ", modlabs[[i]], ": ", ptablelist[[i]]$setModTotals, sep="")),order(paste("Cons ", modlabs$consTOM[[i]], ": ", ptablelist[[i]]$consModTotals, sep=""))],
                     xLabels = paste(" ", modlabs$consTOM[[i]])[order(paste("Cons ", modlabs$consTOM[[i]], ": ", ptablelist[[i]]$consModTotals, sep=""))],
                     yLabels = paste(" ", modlabs[[i]])[order(paste(i, " ", modlabs[[i]], ": ", ptablelist[[i]]$setModTotals, sep=""))],
                     colorLabels = TRUE,
                     xSymbols = paste("Cons ", modlabs$consTOM[[i]], ": ", ptablelist[[i]]$consModTotals, sep="")[order(paste("Cons ", modlabs$consTOM[[i]], ": ", ptablelist[[i]]$consModTotals, sep=""))],
                     ySymbols = paste(i, " ", modlabs[[i]], ": ", ptablelist[[i]]$setModTotals, sep="")[order(paste(i, " ", modlabs[[i]], ": ", ptablelist[[i]]$setModTotals, sep=""))],
                     # textMatrix = ptablelist[[i]]$CountTbl[order(paste(i, " ", modlabs[[i]], ": ", ptablelist[[i]]$setModTotals, sep="")),order(paste("Cons ", modlabs$consTOM[[i]], ": ", ptablelist[[i]]$consModTotals, sep=""))],
                     colors = blueWhiteRed(100)[50:100],
                     main = paste0("Correspondence of ", i, " set-specific and ", i, '-', names(ptablelist)[!names(ptablelist)%in%i], " consensus modules"),
                     cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)
      dev.off()
    })
  }
  
  if (any(runnew['traits'],runnew['selectedTraits'])) {
    cat('Module associations!\n')
    # sampleTrees <- readRDS(file.path(outpath, 'RDS', 'sampleTrees.rds'))
    traits <- readRDS(file.path(outpath, 'RDS', '/traits.rds'))
    multiExpr <- readRDS(file.path(outpath, 'RDS', 'multiExpr_afterscale.rds'))
    matlist <- readRDS(file.path(outpath, 'RDS', 'matlist.rds'))
    nGenes = as.list(rep(checkSets(multiExpr)$nGenes, times = length(matlist)))
    names(nGenes) <- names(matlist)
    nSamples = lapply(matlist, function(i) if(is.integer(nrow(i))) nrow(i) else checkSets(i)$nSamples)
    if (!is.null(mat2)) {
      names(nSamples$consTOM) <- setLabels
    }
    MEslist <- readRDS(file.path(outpath, 'RDS', 'mergedMEs.rds'))
    samp <- as.data.frame.matrix(table(rownames(traits), rownames(traits)))
    traitlist <- list(samp, traits)
    names(traitlist) <- c('Sample', 'Trait')
    orderlists <- vector(mode = 'list', length = length(traitlist))
    names(orderlists) <- names(traitlist)
    # samp[rownames(mat),rownames(mat)]
    for (x in names(traitlist)) {
      cat(x, '\n')
      orderlists[[x]] <- vector(mode = 'list', length = length(MEslist))
      names(orderlists[[x]]) <- names(MEslist)
      for (i in names(MEslist)) {
        cat(i, '\n')
        orderlists[[x]][[i]] <- vector(mode = 'list', length = 2)
        names(orderlists[[x]][[i]]) <- c('traitorder', 'sampleorder')
        if (class(MEslist[[i]])=='list') {
          # cat('List!\n')
          moduleSampCor <- lapply(MEslist[[i]], function(j) cor(j, traitlist[[x]], use = 'p'))
          moduleSampPvalue <- lapply(names(moduleSampCor), function(j) corPvalueStudent(moduleSampCor[[j]], nSamples[[i]][[j]]))
          names(moduleSampPvalue) <- names(moduleSampCor)
          for (j in names(moduleSampCor)) {
            moduleSampCor[[j]][is.na(moduleSampCor[[j]])] <- 0
          }
          if (identical(cor(t(moduleSampCor[[1]])), cor(t(moduleSampCor[[2]])))) {
            modulesampdist <- as.dist(1-cor(t(moduleSampCor[[1]]), use = 'p'))
          } else {
            modulesampdist <- as.dist(1-cor(cor(t(moduleSampCor[[1]])), cor(t(moduleSampCor[[2]])), use = 'p'))
          }
          modulesamphclust <- hclust(modulesampdist, method = agglomeration)
          modulesamphclust <- dendsort(modulesamphclust, type = 'min')
          modulesamporder <- modulesamphclust$labels[modulesamphclust$order]
          orderlists[[x]][[i]][['sampleorder']] <- modulesamporder
          # modulesampdist <- lapply(moduleSampCor, function(j) as.dist(1-cor(t(j), use = 'p')))
          # modulesamphclust <- lapply(modulesampdist, function(j) hclust(j, method = agglomeration))
          # modulesamphclust <- lapply(modulesamphclust, function(j) dendsort(j, type = 'min'))
          # modulesamporder <- lapply(modulesamphclust, function(j) j$labels[j$order])
          # traitsamphclust <- sampleTrees[[i]]
          # traitsamporder <- traitsamphclust$labels[traitsamphclust$order]
          if (identical(cor(t(moduleSampCor[[1]])), cor(t(moduleSampCor[[2]])))) {
            traitsampdist <- as.dist(1-cor(moduleSampCor[[1]], use = 'p'))
          } else {
            traitsampdist <- as.dist(1-cor(cor(moduleSampCor[[1]]), cor(moduleSampCor[[2]]), use = 'p'))
          }
          traitsamphclust <- hclust(traitsampdist, method = agglomeration)
          traitsamphclust <- dendsort(traitsamphclust, type = 'min')
          traitsamporder <- traitsamphclust$labels[traitsamphclust$order]
          orderlists[[x]][[i]][['traitorder']] <- traitsamporder
          for (m in names(moduleSampCor)) {
            moduleSampCor[[m]] <- moduleSampCor[[m]][modulesamporder,traitsamporder]
          }
          for (m in names(moduleSampPvalue)) {
            moduleSampPvalue[[m]] <- moduleSampPvalue[[m]][modulesamporder,traitsamporder]
          }
          # lapply(names(modulesamporder), function(j) saveRDS(modulesamporder[[j]], file.path(outpath, 'RDS', paste0('modulesamporder_', i, '_', j, '.rds')), compress = T))
          # lapply(names(traitsamporder), function(j) saveRDS(traitsamporder[[j]], file.path(outpath, 'RDS', paste0('traitsamporder_', i, '_', j, '.rds')), compress = T))
          
          # Will display correlations and their p-values
          textsampMatrix  <- lapply(names(moduleSampCor), function(j) paste(signif(moduleSampCor[[j]][modulesamporder,traitsamporder], 2), "\n(",
                                                                            signif(moduleSampPvalue[[j]][modulesamporder,traitsamporder], 1), ")", sep = ""))
          names(textsampMatrix) <- names(moduleSampCor)
          for (j in names(textsampMatrix)) {
            dim(textsampMatrix[[j]]) <- dim(moduleSampCor[[j]])
          }
          names(textsampMatrix) <- names(moduleSampCor)
          for (j in names(textsampMatrix)) {
            rownames(textsampMatrix[[j]]) <- modulesamporder
            colnames(textsampMatrix[[j]]) <- traitsamporder
          }
          modcorrsamp <- lapply(names(moduleSampCor), function(j) data.table(melt(moduleSampCor[[j]][modulesamporder,traitsamporder], varnames = c('Module', 'Sample'), value.name = 'Correlation')))
          for (m in modcorrsamp) {
            try(setnames(m, 'value', 'Correlation'), silent = T)
          }
          names(modcorrsamp) <- names(moduleSampCor)
          modpvalsamp <- lapply(names(moduleSampPvalue), function(j) data.table(melt(moduleSampPvalue[[j]][modulesamporder,traitsamporder], varnames = c('Module', 'Sample'), value.name = 'p value')))
          for (m in modpvalsamp) {
            try(setnames(m, 'value', 'p value'), silent = T)
          }
          names(modpvalsamp) <- names(moduleSampPvalue)
          modcorexportsamp <- lapply(names(modcorrsamp), function(j) modcorexportsamp <- merge(modcorrsamp[[j]], modpvalsamp[[j]], by=c('Module', 'Sample'), all = T))
          names(modcorexportsamp) <- names(modcorrsamp)
          lapply(modcorexportsamp, function(j) {
            j[,Correlation:=signif(Correlation, 3)]
            j[,`p value`:=signif(`p value`, 3)]
          })
          lapply(names(modcorexportsamp), function(j) fwrite(modcorexportsamp[[j]], file.path(outpath, 'CSVs', paste0('Module_', x, '_Correlation_',i, '_', j, '.csv')), sep = ',', col.names = T, row.names = F))
          # Display the correlation values within a heatmap plot
          lapply(names(moduleSampCor), function(j) {
            png(file.path(outpath, 'PLOTS', paste0('module_', x, '_correlations_', i, '_', j, '.png')), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
            par(mar = c(6, 5, 3, 1))
            labeledHeatmap(Matrix = moduleSampCor[[j]][modulesamporder,traitsamporder], cex.lab.y = 0.6, cex.lab.x = 0.4, yColorWidth = 0.01,
                           xLabels = traitsamporder,
                           yLabels = modulesamporder,
                           ySymbols = modulesamporder,
                           colorLabels = FALSE,
                           colors = colorRampPalette(c(TUMblue, 'white', TUMred))(50),
                           # textMatrix = textMatrix,
                           setStdMargins = FALSE,
                           cex.text = 0.35,
                           zlim = c(-1,1),
                           main = paste0("Module-", x, " relationships\n", i, ' - ', j))
            dev.off()
          })
          # HERE HERE HERE HERE
          consensusCor = matrix(NA, nrow(moduleSampCor[[1]]), ncol(moduleSampCor[[1]]))
          rownames(consensusCor) <- rownames(moduleSampCor[[1]])
          colnames(consensusCor) <- colnames(moduleSampCor[[1]])
          consensusPvalue = matrix(NA, nrow(moduleSampCor[[1]]), ncol(moduleSampCor[[1]]))
          rownames(consensusPvalue) <- rownames(moduleSampCor[[1]])
          colnames(consensusPvalue) <- colnames(moduleSampCor[[1]])
          # Find consensus negative correlations
          negative = moduleSampCor[[1]] < 0 & moduleSampCor[[2]] < 0
          consensusCor[negative] = pmax(moduleSampCor[[1]][negative], moduleSampCor[[2]][negative])
          consensusPvalue[negative] = pmax(moduleSampPvalue[[1]][negative], moduleSampPvalue[[2]][negative])
          # Find consensus positive correlations
          positive = moduleSampCor[[1]] > 0 & moduleSampCor[[2]] > 0
          consensusCor[positive] = pmin(moduleSampCor[[1]][positive], moduleSampCor[[2]][positive])
          consensusPvalue[positive] = pmax(moduleSampPvalue[[1]][positive], moduleSampPvalue[[2]][positive])
          textMatrix = paste(signif(consensusCor[modulesamporder,traitsamporder], 2), "\n(",
                             signif(consensusPvalue[modulesamporder,traitsamporder], 1), ")", sep = "");
          dim(textMatrix) = dim(consensusCor)
          #pdf(file = "Plots/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
          png(file.path(outpath, 'PLOTS', paste0('module_', x, '_correlations_', i, '.png')), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
          par(mar = c(6, 5, 3, 1))
          labeledHeatmap(Matrix = consensusCor[modulesamporder,traitsamporder], cex.lab.y = 0.6, cex.lab.x = 0.4, yColorWidth = 0.01,
                         xLabels = traitsamporder,
                         yLabels = modulesamporder,
                         ySymbols = modulesamporder,
                         colorLabels = FALSE,
                         colors = colorRampPalette(c(TUMblue, 'white', TUMred))(50),
                         textMatrix = NULL, # textMatrix
                         setStdMargins = FALSE,
                         cex.text = 0.35,
                         zlim = c(-1,1),
                         main = paste("Consensus module-", x, " relationships across\n",
                                      paste(setLabels, collapse = " and "), sep = ''))
          dev.off()
        } else {
          # cat('Vec!\n')
          moduleSampCor <- cor(MEslist[[i]], traitlist[[x]], use = 'p')
          moduleSampPvalue = corPvalueStudent(moduleSampCor, nSamples[[i]])
          modulesampdist <- as.dist(1-cor(t(moduleSampCor), use = 'p'))
          modulesamphclust <- hclust(modulesampdist, method = agglomeration)
          modulesamphclust <- dendsort(modulesamphclust, type = 'min')
          modulesamporder <- modulesamphclust$labels[modulesamphclust$order]
          orderlists[[x]][[i]][['sampleorder']] <- modulesamporder
          traitsampdist <- as.dist(1-cor(moduleSampCor, use = 'p'))
          traitsamphclust <- hclust(traitsampdist, method = agglomeration)
          traitsamphclust <- dendsort(traitsamphclust, type = 'min')
          traitsamporder <- traitsamphclust$labels[traitsamphclust$order]
          orderlists[[x]][[i]][['traitorder']] <- traitsamporder
          moduleSampCor <- moduleSampCor[modulesamporder,traitsamporder]
          moduleSampPvalue <- moduleSampPvalue[modulesamporder,traitsamporder]
          # saveRDS(modulesamporder, file.path(outpath, 'RDS', paste0('modulesamporder_', i, '.rds')), compress = T)
          # saveRDS(traitsamporder, file.path(outpath, 'RDS', paste0('traitsamporder_', i, '.rds')), compress = T)
          # Will display correlations and their p-values
          textsampMatrix =  paste(signif(moduleSampCor[modulesamporder,traitsamporder], 2), "\n(",
                                  signif(moduleSampPvalue[modulesamporder,traitsamporder], 1), ")", sep = "")
          dim(textsampMatrix) = dim(moduleSampCor)
          rownames(textsampMatrix) <- modulesamporder
          colnames(textsampMatrix) <- traitsamporder
          modcorrsamp <- data.table(melt(moduleSampCor[modulesamporder,traitsamporder], varnames = c('Module', 'Sample'), value.name = 'Correlation'))
          try(setnames(modcorrsamp, 'value', 'Correlation'), silent = T)
          modpvalsamp <- data.table(melt(moduleSampPvalue[modulesamporder,traitsamporder], varnames = c('Module', 'Sample'), value.name = 'p value'))
          try(setnames(modpvalsamp, 'value', 'p value'), silent = T)
          modcorexportsamp <- merge(modcorrsamp, modpvalsamp, by=c('Module', 'Sample'), all = T)
          modcorexportsamp[,Correlation:=signif(Correlation, 3)]
          modcorexportsamp[,`p value`:=signif(`p value`, 3)]
          fwrite(modcorexportsamp, file.path(outpath, 'CSVs', paste0('Module_', x, '_Correlation_', i, '.csv')), sep = ',', col.names = T, row.names = F)
          # Display the correlation values within a heatmap plot
          png(file.path(outpath, 'PLOTS', paste0('module_', x, '_correlations_', i, '.png')), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
          par(mar = c(6, 5, 3, 1))
          labeledHeatmap(Matrix = moduleSampCor[modulesamporder,traitsamporder], cex.lab.y = 0.6, cex.lab.x = 0.4, yColorWidth = 0.01,
                         xLabels = traitsamporder,
                         yLabels = modulesamporder,
                         ySymbols = modulesamporder,
                         colorLabels = FALSE,
                         colors = colorRampPalette(c(TUMblue, 'white', TUMred))(50),
                         # textMatrix = textMatrix,
                         setStdMargins = FALSE,
                         cex.text = 0.35,
                         zlim = c(-1,1),
                         main = paste0("Module-", x, " relationships\n", i))
          dev.off()
        }
      }
      # return(list(traitorderlist=traitorderlist, sampleorderlist=sampleorderlist))
      saveRDS(orderlists, file.path(outpath, 'RDS', 'orderlists.rds'), compress = T)
    }
      
      
      
    
    # tmp <- traits[colnames(mat1),selectedTraits,drop=F]
    # tmp <- tmp[,colSums(tmp, na.rm = T)!=0]
    # moduleTraitCor = cor(MEs, tmp, use = "p") # stats:::    %in%names(which(colSums(traits[rownames(mat),selectedTraits,drop=F])!=0))
    # moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    # 
    # moduledist <- as.dist(1-cor(t(moduleTraitCor), use = 'p'))
    # modulehclust <- hclust(moduledist, method = agglomeration)
    # modulehclust <- dendsort(modulehclust, type = 'min')
    # moduleorder <- modulehclust$labels[modulehclust$order]
    # traitdist <- as.dist(1-cor(moduleTraitCor, use = 'p'))
    # traithclust <- hclust(traitdist, method = agglomeration)
    # traithclust <- dendsort(traithclust, type = 'min')
    # traitorder <- traithclust$labels[traithclust$order]
    # 
    # saveRDS(moduleorder, file.path(outpath, 'RDS', 'moduleorder.rds'), compress = T)
    # saveRDS(traitorder, file.path(outpath, 'RDS', 'traitorder.rds'), compress = T)
    # 
    # # Will display correlations and their p-values
    # textMatrix =  paste(signif(moduleTraitCor[moduleorder,traitorder], 2), "\n(",
    #                     signif(moduleTraitPvalue[moduleorder,traitorder], 1), ")", sep = "")
    # dim(textMatrix) = dim(moduleTraitCor)
    # rownames(textMatrix) <- moduleorder
    # colnames(textMatrix) <- traitorder
    # modcorr <- data.table(melt(moduleTraitCor[moduleorder,traitorder], varnames = c('Module', 'Trait'), value.name = 'Correlation'))
    # modpval <- data.table(melt(moduleTraitPvalue[moduleorder,traitorder], varnames = c('Module', 'Trait'), value.name = 'p value'))
    # modcorexport <- merge(modcorr, modpval, by=c('Module', 'Trait'), all = T)
    # modcorexport[,Correlation:=signif(Correlation, 3)]
    # modcorexport[,`p value`:=signif(`p value`, 3)]
    # fwrite(modcorexport, file.path(outpath, 'CSVs', 'ModuleCorrelation.csv'), sep = ',', col.names = T, row.names = F)
    # 
    # # Display the correlation values within a heatmap plot
    # png(file.path(outpath, 'PLOTS', 'module_associations.png'), width = (210-32.5*2), height = (210-32.5*2), units = "mm", res = 300, pointsize = 10)
    # par(mar = c(6, 5, 3, 1))
    # labeledHeatmap(Matrix = moduleTraitCor[moduleorder,traitorder], cex.lab.y = 0.6, cex.lab.x = 0.4, yColorWidth = 0.01,
    #                xLabels = traitorder,
    #                yLabels = moduleorder,
    #                ySymbols = moduleorder,
    #                colorLabels = FALSE,
    #                colors = colorRampPalette(c(TUMblue, 'white', TUMred))(50),
    #                # textMatrix = textMatrix,
    #                setStdMargins = FALSE,
    #                cex.text = 0.35,
    #                zlim = c(-1,1),
    #                main = paste("Module-trait relationships"))
    # dev.off()  
  }
  
  # traitorder <- readRDS(file.path(outpath, 'RDS', 'traitorder.rds'))
  matlist <- readRDS(file.path(outpath, 'RDS', 'matlist.rds'))
  orderlists <- readRDS(file.path(outpath, 'RDS', 'orderlists.rds'))
  geneTreelist <- readRDS(file.path(outpath, 'RDS', 'geneTree.rds'))
  moduleColors <- readRDS(file.path(outpath, 'RDS', 'mergedColors.rds'))
  mergedColors <- moduleColors
  anndtlist <-lapply(names(geneTreelist), function(i) data.table(rn=geneTreelist[[i]]$labels[geneTreelist[[i]]$order], ord=seq(1:length(mergedColors[[i]])), col=mergedColors[[i]][geneTreelist[[i]]$order]))
  names(anndtlist) <- names(geneTreelist)
  
  options(warn = -1)
  if (is.null(colord)) {
    potcols <- sort(unique(mergedColors))
    colord <- c()
    while (length(potcols[!potcols%in%colord])>0) {
      potcols <- potcols[!potcols%in%colord]
      idx <- as.numeric(readline(prompt = paste0('Select first colour of first module:\n', paste0(paste0(potcols, ' [', seq_along(potcols), ']', collapse = '\n'), '\n'))))
      if (!is.na(idx) && idx%%1==0 && idx<=length(potcols)) {
        colord[length(colord)+1] <- potcols[idx]
      } else {
        cat('Please select a number between 1 and', length(potcols))
      }
    }
  }
  options(warn = 0)
  
  if (colord == 'hclust') {
    # MEorder <- readRDS(file.path(outpath, 'RDS', 'MEorder.rds'))
    # colord <- gsub('^ME', '', MEorder)
    # moduleorder <- readRDS(file.path(outpath, 'RDS', 'moduleorder.rds'))
    colordlist <- lapply(orderlists$Trait, function(i) gsub('^ME', '', i$sampleorder))
  }
  
  coltoplotlist <- NULL
  if (colord == 'sorted') {
    # moduleorder <- readRDS(file.path(outpath, 'RDS', 'moduleorder.rds'))
    # orderlist <- mclapply(
    #   gsub(
    #     '^ME',
    #     '',
    #     moduleorder),
    #   function(i)
    #     rownames(mat)[order(
    #       unlist(
    #         lapply(
    #           rownames(mat),
    #           function(j)
    #             sum(
    #               rowMedians(
    #                 mat[
    #                   rownames(mat)==j,
    #                   anndt[col==i,rn],
    #                   drop=F],
    #                 na.rm = T)
    #               /
    #                 rowMedians(
    #                   mat[
    #                     rownames(mat)!=j,
    #                     anndt[col==i,rn],
    #                     drop = F], na.rm = T)))), decreasing = T)],
    #   mc.cores = 10)
    # names(orderlist) <- gsub(
    #   '^ME',
    #   '',
    #   moduleorder)
    
    orderdtlist <- lapply(names(matlist), function(x) {
      if (class(matlist[[x]])=='list') {
        # cat('here\n')
        lapply(matlist[[x]], function(m) {
          do.call(rbind, mclapply(rownames(m$data), function(i)
            sapply(lapply(orderlists$Trait, function(u) gsub('^ME', '', u$sampleorder))[[x]], function(j)
              median(unlist(m$data[i,anndtlist[[x]][col==j,rn]]), na.rm = T)), mc.cores = 10))
        })
      } else {
        do.call(rbind, mclapply(rownames(matlist[[x]]), function(i)
        sapply(lapply(orderlists$Trait, function(u) gsub('^ME', '', u$sampleorder))[[x]], function(j)
          median(unlist(matlist[[x]][i,anndtlist[[x]][col==j,rn]]), na.rm = T)), mc.cores = 10))
      }
    })
    
    names(orderdtlist) <- names(matlist)
    
    # orderdt <- do.call(rbind, mclapply(rownames(mat), function(i)
    #   sapply(gsub('^ME', '', moduleorder), function(j)
    #     median(mat[i,anndt[col==j,rn]], na.rm = T)), mc.cores = 10))
    
    for (i in names(orderdtlist)) {
      if (class(orderdtlist[[i]])=='list') {
        for (j in names(orderdtlist[[i]])) {
          rownames(orderdtlist[[i]][[j]]) <- rownames(matlist[[i]][[j]]$data)
          orderdtlist[[i]][[j]] <- t(scale(t(orderdtlist[[i]][[j]])))    
        }  
      } else {
        rownames(orderdtlist[[i]]) <- rownames(matlist[[i]])
        # saveRDS(orderdtlist, 'sec/FL/orderlistdebug.rds', compress = T)
        orderdtlist[[i]] <- t(scale(t(orderdtlist[[i]])))
      }
    }
    
    toassignlist <- lapply(orderdtlist, function(i) {
      if (class(i)=='list') {
        lapply(i, function(j) colnames(j))  
      } else {
        colnames(i)
      }
      })
    usemax <- T
    mxlist <- lapply(orderdtlist, function(i) {
      if(class(i)=='list') {
        lapply(i, function(j) min(apply(j, 1, max, na.rm = T)))
      } else {
        min(apply(i, 1, max, na.rm = T))
      }
    })
    # mx <- min(apply(orderdt, 1, function(i) max(i, na.rm = T)))
    reslist <- vector(mode = 'list', length = length(names(toassignlist)))
    names(reslist) <- names(toassignlist)
    for (m in names(toassignlist)) {
      cat(m, '\n')
      if (class(toassignlist[[m]])=='list') {
        reslist[[m]] <- vector(mode = 'list', length = length(names(toassignlist[[m]])))
        names(reslist[[m]]) <- names(toassignlist[[m]])
        for (k in names(toassignlist[[m]])) {
          cat(k, '\n')
          reslist[[m]][[k]] <- list()
          usemax <- T
          while (length(toassignlist[[m]][[k]])>0) {
            if (usemax) {
              for (i in rownames(orderdtlist[[m]][[k]])) {
                reslist[[m]][[k]][[i]] <- append(reslist[[m]][[k]][[i]], if(max(orderdtlist[[m]][[k]][i,colnames(orderdtlist[[m]][[k]])%in%toassignlist[[m]][[k]]], na.rm = T) >= mxlist[[m]][[k]]) colnames(orderdtlist[[m]][[k]][i,colnames(orderdtlist[[m]][[k]])%in%toassignlist[[m]][[k]],drop=F][1,which.max(orderdtlist[[m]][[k]][i,colnames(orderdtlist[[m]][[k]])%in%toassignlist[[m]][[k]],drop=F]),drop=F]) else NA)
              }
            } else {
              for (i in rownames(orderdtlist[[m]][[k]])) {
                reslist[[m]][[k]][[i]] <- append(reslist[[m]][[k]][[i]], colnames(orderdtlist[[m]][[k]][i,colnames(orderdtlist[[m]][[k]])%in%toassignlist[[m]][[k]],drop=F][1,which.max(orderdtlist[[m]][[k]][i,colnames(orderdtlist[[m]][[k]])%in%toassignlist[[m]][[k]],drop=F]),drop=F]))
              }
            }
            assigned <- unique(unlist(reslist[[m]][[k]]))
            toassignlist[[m]][[k]] <- toassignlist[[m]][[k]][!toassignlist[[m]][[k]]%in%assigned]
            if((length(unique(unlist(lapply(reslist[[m]][[k]], function(i) tail(i, 1)))))==1 && is.na(unique(unlist(lapply(reslist[[m]][[k]], function(i) tail(i, 1)))))) & length(toassignlist[[m]][[k]])!=0) {
              usemax <- F
            }
            cat(usemax, '\n')
          }
          reslist[[m]][[k]] <- do.call(rbind, reslist[[m]][[k]])
        }
      } else {
        reslist[[m]] <- list()
        usemax <- T
        while (length(toassignlist[[m]])>0) {
          if (usemax) {
            for (i in rownames(orderdtlist[[m]])) {
              reslist[[m]][[i]] <- append(reslist[[m]][[i]], if(max(orderdtlist[[m]][i,colnames(orderdtlist[[m]])%in%toassignlist[[m]]], na.rm = T) >= mxlist[[m]]) colnames(orderdtlist[[m]][i,colnames(orderdtlist[[m]])%in%toassignlist[[m]],drop=F][1,which.max(orderdtlist[[m]][i,colnames(orderdtlist[[m]])%in%toassignlist[[m]],drop=F]),drop=F]) else NA)
            }
          } else {
            for (i in rownames(orderdtlist[[m]])) {
              reslist[[m]][[i]] <- append(reslist[[m]][[i]], colnames(orderdtlist[[m]][i,colnames(orderdtlist[[m]])%in%toassignlist[[m]],drop=F][1,which.max(orderdtlist[[m]][i,colnames(orderdtlist[[m]])%in%toassignlist[[m]],drop=F]),drop=F]))
            }
          }
          assigned <- unique(unlist(reslist[[m]]))
          toassignlist[[m]] <- toassignlist[[m]][!toassignlist[[m]]%in%assigned]
          if((length(unique(unlist(lapply(reslist[[m]], function(i) tail(i, 1)))))==1 && is.na(unique(unlist(lapply(reslist[[m]], function(i) tail(i, 1)))))) & length(toassignlist[[m]])!=0) {
            usemax <- F
          }
          cat(usemax, '\n')
        }
        reslist[[m]] <- do.call(rbind, reslist[[m]])
      }
    }
    
    colordlist <- lapply(reslist, function(i) {
      if (class(i)=='list') {
        lapply(i, function(j) unique(as.character(melt(j)$value)))
      } else {
        unique(as.character(melt(i)$value))  
      }
    })
    for (i in names(colordlist)) {
      if (class(colordlist[[i]])=='list') {
        for (j in names(colordlist[[i]])) {
          colordlist[[i]][[j]] <- colordlist[[i]][[j]][!is.na(colordlist[[i]][[j]])]
          colordlist[[i]][[j]] <- rev(colordlist[[i]][[j]])
        }
      } else {
        colordlist[[i]] <- colordlist[[i]][!is.na(colordlist[[i]])]
        colordlist[[i]] <- rev(colordlist[[i]])
      }
    }
    # colord <- colord[!is.na(colord)]
    # colord <- rev(colord)
    seldtlist <- lapply(reslist, function(i) {
      if (class(i)=='list') {
        lapply(i, function(j) data.table(melt(j)))
      } else {
        data.table(melt(i))
      }
    })
    for (i in seldtlist) {
      if (all(class(i)=='list')) {
        for (j in i) {
          try(setnames(j, c('X1', 'X2'), c('Var1', 'Var2')))
        }
      } else {
        try(setnames(i, c('X1', 'X2'), c('Var1', 'Var2')))
      }
    }
    # seldt <- data.table(melt(do.call(rbind, reslist)))
    # rm(reslist)
    # gc()
    potphaselist <- lapply(seldtlist, function(i) {
      if (any(class(i)=='list')) {
        lapply(i, function(j) j[,length(unique(value)),by=Var2][V1%in%1,Var2])
      } else {
        i[,length(unique(value)),by=Var2][V1%in%1,Var2]
      }
    })
    # potphase <- seldt[,length(unique(value)),by=Var2][V1%in%1,Var2]
    coltoplotlist <- vector(mode = 'list', length = length(names(seldtlist)))
    names(coltoplotlist) <- names(seldtlist)
    for (i in names(seldtlist)) {
      if (any(class(seldtlist[[i]])=='list')) {
        coltoplotlist[[i]] <- vector(mode = 'list', length = length(names(seldtlist[[i]])))
        names(coltoplotlist[[i]]) <- names(seldtlist[[i]])
        for (j in names(seldtlist[[i]])) {
          coltoplotlist[[i]][[j]] <- seldtlist[[i]][[j]][Var2<=seldtlist[[i]][[j]][Var2%in%potphaselist[[i]][[j]],unique(value),by=Var2][is.na(V1),Var2-1],unique(value)]
        }
      } else {
        coltoplotlist[[i]] <- seldtlist[[i]][Var2<=seldtlist[[i]][Var2%in%potphaselist[[i]],unique(value),by=Var2][is.na(V1),Var2-1],unique(value)]
      }
    }
    # coltoplot <- seldt[Var2<=seldt[Var2%in%potphase,unique(value),by=Var2][is.na(V1),Var2-1],unique(value)]
    for (i in names(coltoplotlist)) {
      if (class(coltoplotlist[[i]])=='list') {
        for (j in names(coltoplotlist[[i]])) {
          coltoplotlist[[i]][[j]] <- coltoplotlist[[i]][[j]][!is.na(coltoplotlist[[i]][[j]])]
        }
      } else {
        coltoplotlist[[i]] <- coltoplotlist[[i]][!is.na(coltoplotlist[[i]])]
      }
    }
    # coltoplot <- coltoplot[!is.na(coltoplot)]
    # coltoplot <- as.character(coltoplot)
    for (i in names(potphaselist)) {
      if (class(potphaselist[[i]])=='list') {
        for (j in names(potphaselist[[i]])) {
          if (length(potphaselist[[i]][[j]])==0) {
            coltoplotlist[[i]][[j]] <- gsub('^ME', '', orderlists$Trait[[i]]$sampleorder)
          }
        }
      } else {
        if (length(potphaselist[[i]])==0) {
          coltoplotlist[[i]] <- gsub('^ME', '', orderlists$Trait[[i]]$sampleorder)
        }
      }
    }
    # invisible(ifelse(length(potphase)==0,coltoplot <- gsub('^ME', '', moduleorder),coltoplot <- coltoplot))
    
    
    # tt <- do.call(rbind, lapply(orderlist, function(i) match(roword, i)))
    # colnames(tt) <- roword
    # # colord <- names(sort(rowSums(scale(tt))))
    # windowmin <- function(win=4, dt, slide = T) {
    #   if (slide) {
    #     winlist <- lapply(1:(ncol(dt)-(win-1)), function(i) seq(i, i+(win-1)))
    #   } else {
    #     ct <- as.character(cut(1:ncol(dt), seq(0, ncol(dt), win)))
    #     ct[is.na(ct)] <- 'last'
    #     winlist <- lapply(unique(ct), function(i) seq(1:ncol(dt))[ct%in%i])
    #   }
    #   nm <- do.call(cbind, lapply(winlist, function(i) rowMax(dt[,i,drop=F])))
    #   rownames(nm) <- rownames(dt)
    #   # tmp <- apply(nm, 2, rescale , c(0,1))
    #   # tmp <- nm
    #   tmp <- scale(nm)
    #   rownames(tmp) <- rownames(nm)
    #   tmp <-  rowSums(tmp)
    #   colord <- names(sort(tmp))
    #   return(colord)
    # }
    # colord <- windowmin(win = 3, dt = tt, slide = T)
  }
  
  colannlist <- vector(mode = 'list', length = length(names(matlist)))
  names(colannlist) <- names(matlist)  
  for (i in names(matlist)) {
    if (class(matlist[[i]])=='list') {
      colannlist[[i]] <- vector(mode = 'list', length = length(names(matlist[[i]])))
      names(colannlist[[i]]) <- names(matlist[[i]])
      for (j in names(matlist[[i]])) {
        colannlist[[i]][[j]] <- data.frame(row.names = colnames(matlist[[i]][[j]]$data), ModuleMembership = moduleColors[[i]])
      }
    } else {
      colannlist[[i]] <- data.frame(row.names = colnames(matlist[[i]]), ModuleMembership = moduleColors[[i]]) 
    }
  }
  # colann <- data.frame(row.names = colnames(mat), ModuleMembership=moduleColors, stringsAsFactors = F)
  ModuleMembershiplist <- lapply(colannlist, function(i) {
    if (class(i)=='list') {
      unique(i[[1]]$ModuleMembership)
    } else {
      unique(i$ModuleMembership)
    }
  })
  for (i in names(ModuleMembershiplist)) {
    # if (class(ModuleMembershiplist[[i]])=='list') {
    #   for (j in names(ModuleMembershiplist[[i]])) {
    #     names(ModuleMembershiplist[[i]][[j]]) <- ModuleMembershiplist[[i]][[j]]
    #   }
    # } else {
      names(ModuleMembershiplist[[i]]) <- ModuleMembershiplist[[i]]
    # }
  }
  # ModuleMembership <- unique(colann$ModuleMembership)
  # names(ModuleMembership) <- unique(colann$ModuleMembership)
  
  if (heatmaps) {
    cat('Heatmaps!\n')
    traits <- readRDS(file.path(outpath, 'RDS', 'traits.rds'))
    matlist <- readRDS(file.path(outpath, 'RDS', 'matlist.rds'))
    if (length(sampleord)==1 && sampleord=='hclust') {
      sampleTrees <- readRDS(file.path(outpath, 'RDS', 'sampleTrees.rds'))
      # sampleTree2 <- readRDS(file.path(outpath, 'RDS', 'sampleTree.rds'))
      rowordlist <- lapply(sampleTrees, function(i) i$labels[i$order])
      # roword <- sampleTree2$labels[sampleTree2$order]
      for (i in names(matlist)) {
        if (class(matlist[[i]])=='list') {
          for (j in names(matlist[[i]])) {
            matlist[[i]][[j]]$data <- matlist[[i]][[j]]$data[rowordlist[[i]],]
          }
        } else {
          matlist[[i]] <- matlist[[i]][rowordlist[[i]],]
        }
      }
      # mat <- mat[roword,]
    } else {
      sampleTrees <- readRDS(file.path(outpath, 'RDS', 'sampleTrees.rds'))
      rowordlist <- lapply(sampleTrees, function(i) sampleord)
      for (i in names(matlist)) {
        if (class(matlist[[i]])=='list') {
          for (j in names(matlist[[i]])) {
            if (all(!is.na(match(rownames(matlist[[i]][[j]]$data, sampleord))))) {
              roword <- sampleord
              matlist[[i]][[j]]$data <- matlist[[i]][[j]]$data[roword,]
            } else {
              stop('sampleord does not match rownames of matrix!')
            }
          }
        } else {
          if (all(!is.na(match(rownames(matlist[[i]]), sampleord)))) {
            roword <- sampleord
            matlist[[i]] <- matlist[[i]][roword,]
          } else {
            stop('sampleord does not match rownames of matrix!')
          }
        }
      }
      # if (all(!is.na(match(rownames(mat), sampleord)))) {
      #   roword <- sampleord
      #   mat <- mat[roword,]
      # } else {
      #   stop('sampleord does not match rownames of matrix!')
      # }
    }
    
    
    # coldist <- as.dist(1-cor(mat, use = 'p'))
    # rowdist <- as.dist(1-cor(t(mat), use = 'p'))
    rowannlist <- lapply(matlist, function(i) {
      if(class(i)=='list') {
        lapply(i, function(j) {
          tmp <- data.frame(row.names = rownames(j$data), traits[rownames(j$data),plottedTraits,drop=F], check.names = F)
          colnames(tmp) <- plottedTraits
          tmp <- tmp[,colSums(tmp, na.rm = T)!=0,drop = F]
          return(tmp)
        })
      } else {
        tmp <- data.frame(row.names = rownames(i), traits[rownames(i),plottedTraits,drop=F], check.names = F)
        colnames(tmp) <- plottedTraits
        tmp <- tmp[,colSums(tmp, na.rm = T)!=0,drop = F]
        return(tmp)
      }
    })
    
    # rowann <- data.frame(row.names = rownames(mat), traits[rownames(mat),plottedTraits], check.names = F)
    # colnames(rowann) <- plottedTraits

    anncolslist <- lapply(ModuleMembershiplist, function(i) {
      list(ModuleMembership=i)
    })
    # anncols <- list(ModuleMembership=ModuleMembership)
    cat('Here!\n')
    for (i in names(anncolslist)) {
      tmp <- c(TUMturk, TUMbeige)
      names(tmp) <- c('1', '0')
      for (j in colnames(traits)) {
        anncolslist[[i]][[j]] <- tmp
      }  
    }
    cat('Not Here!\n')
    # tmp <- c(TUMturk, TUMbeige)
    # names(tmp) <- c('1', '0')
    # for (i in colnames(traits)) {
    #   anncols[[i]] <- tmp
    # }
    # rowann <- rowann[,colSums(rowann, na.rm = T)!=0,drop = F]
    # rowann[is.na(rowann)] <- 'NA'
    # rowann[rowann=='NA'] <- NA
    
    for (i in anndtlist) {
      if (all(class(i)=='list')) {
        for (j in i) {
          setkey(j, ord)
        }
      } else {
        setkey(i, ord)
      }
    }
    # setkey(anndt, ord)
    # require(pheatmap)
    
    for (i in names(matlist)) {
      if (class(matlist[[i]])=='list') {
        for (j in names(matlist[[i]])) {
          png(file.path(outpath, 'PLOTS', paste0('heatmap_wgcnaorder_', i, '_', j, '.png')), width = (210-32.5*2)*2, height = (210-32.5*2)*2, units = "mm", res = 300, pointsize = 10)
          pheatmap(matlist[[i]][[j]]$data[rowordlist[[i]], # [order(traits$LeafClassNumeric),]
                       anndtlist[[i]]$rn],
                   na_col = 'black',
                   scale = "column",
                   # clustering_distance_rows = rowdist,
                   # clustering_distance_cols = coldist,
                   cluster_rows = F,
                   cluster_cols = F,
                   clustering_method = agglomeration,
                   show_colnames = F,
                   # breaks = seq(
                   #   min(mat[rownames(mat)%in%rownames(traits), # [order(traits$LeafClassNumeric),]
                   #           anndt$rn], na.rm =T),
                   #   max(mat[rownames(mat)%in%rownames(traits), # [order(traits$LeafClassNumeric),]
                   #           anndt$rn], na.rm = T), length.out = 256),
                   color = colorRampPalette(c('black', 'black', 'yellow'))(255),
                   annotation_row = rowannlist[[i]][[j]][,orderlists$Trait[[i]]$traitorder[orderlists$Trait[[i]]$traitorder%in%plottedTraits],drop=F],
                   annotation_col = colannlist[[i]][[j]][,,drop=F],
                   annotation_colors = anncolslist[[i]])
          dev.off()
        }
      } else {
        png(file.path(outpath, 'PLOTS', paste0('heatmap_wgcnaorder_', i, '.png')), width = (210-32.5*2)*2, height = (210-32.5*2)*2, units = "mm", res = 300, pointsize = 10)
        pheatmap(matlist[[i]][rowordlist[[i]], # [order(traits$LeafClassNumeric),]
                                        anndtlist[[i]]$rn],
                 na_col = 'black',
                 scale = "column",
                 # clustering_distance_rows = rowdist,
                 # clustering_distance_cols = coldist,
                 cluster_rows = F,
                 cluster_cols = F,
                 clustering_method = agglomeration,
                 show_colnames = F,
                 # breaks = seq(
                 #   min(mat[rownames(mat)%in%rownames(traits), # [order(traits$LeafClassNumeric),]
                 #           anndt$rn], na.rm =T),
                 #   max(mat[rownames(mat)%in%rownames(traits), # [order(traits$LeafClassNumeric),]
                 #           anndt$rn], na.rm = T), length.out = 256),
                 color = colorRampPalette(c('black', 'black', 'yellow'))(255),
                 annotation_row = rowannlist[[i]][,orderlists$Trait[[i]]$traitorder[orderlists$Trait[[i]]$traitorder%in%plottedTraits],drop=F],
                 annotation_col = colannlist[[i]][,,drop=F],
                 annotation_colors = anncolslist[[i]])
        dev.off()
      }
    }
    # png(file.path(outpath, 'PLOTS', 'heatmap_wgcnaorder.png'), width = (210-32.5*2)*2, height = (210-32.5*2)*2, units = "mm", res = 300, pointsize = 10)
    # pheatmap(mat[roword, # [order(traits$LeafClassNumeric),]
    #              anndt$rn],
    #          na_col = 'black',
    #          # clustering_distance_rows = rowdist,
    #          # clustering_distance_cols = coldist,
    #          cluster_rows = F,
    #          cluster_cols = F,
    #          clustering_method = agglomeration,
    #          show_colnames = F,
    #          # breaks = seq(
    #          #   min(mat[rownames(mat)%in%rownames(traits), # [order(traits$LeafClassNumeric),]
    #          #           anndt$rn], na.rm =T),
    #          #   max(mat[rownames(mat)%in%rownames(traits), # [order(traits$LeafClassNumeric),]
    #          #           anndt$rn], na.rm = T), length.out = 256),
    #          color = colorRampPalette(c('black', 'black', 'yellow'))(255),
    #          annotation_row = rowann[,traitorder[traitorder%in%plottedTraits],drop=F],
    #          annotation_col = colann[,,drop=F],
    #          annotation_colors = anncols)
    # dev.off()
    
    
    if (is.null(coltoplotlist)) {
      coltoplotlist <- lapply(orderlists$Trait, function(i) gsub('^ME', '', i$sampleorder))
      # coltoplot <- gsub('^ME', '', moduleorder)
    }
    
    if (!is.null(mat2)) {
      tmp <- anndtlist[["consTOM"]]
      anndtlist[["consTOM"]] <- vector(mode = 'list', length = length(setLabels))
      names(anndtlist[['consTOM']]) <- setLabels
      for (i in setLabels) {
        anndtlist[['consTOM']][[i]] <- tmp
      }
    }
    for (i in names(anndtlist)) {
      if (all(class(anndtlist[[i]])=='list')) {
        for (j in names(anndtlist[[i]])) {
          anndtlist[[i]][[j]][,col:=factor(col, levels = colordlist[[i]][[j]], ordered = T, labels = colordlist[[i]][[j]])]
          setkey(anndtlist[[i]][[j]], col, ord)
        }
      } else {
        anndtlist[[i]][,col:=factor(col, levels = colordlist[[i]], ordered = T, labels = colordlist[[i]])]
        setkey(anndtlist[[i]], col, ord)
      }
    }
    # anndt[,col:=factor(col, levels = colord, ordered = T, labels = colord)]
    # setkey(anndt, col, ord)
    
    # scalecut <- 0
    for (i in seq_along(matlist)) {
      if (all(class(matlist[[i]])=='list')) {
        for (j in seq_along(matlist[[i]])) {
          matlist[[i]][[j]]$data <- scale(matlist[[i]][[j]]$data)
        }
      } else {
        matlist[[i]] <- scale(matlist[[i]])
      }
    }
    lowerlist <- lapply(matlist, function(i) {
      if (class(i)=='list') {
        lapply(i, function(j) {
          as.numeric(unlist(lapply(strsplit(levels(cut(j$data[j$data<scalecut], breaks = 2, na.rm = T, dig.lab = -1)), '\\(|\\,|\\]'), function(m) m[[2]])))
        })
      } else {
        as.numeric(unlist(lapply(strsplit(levels(cut(i[i<scalecut], breaks = 2, na.rm = T, dig.lab = -1)), '\\(|\\,|\\]'), function(m) m[[2]])))
      }
    })
    higherlist <- lapply(matlist, function(i) {
      if (class(i)=='list') {
        lapply(i, function(j) {
          as.numeric(unlist(lapply(strsplit(levels(cut(j$data[j$data>=scalecut], breaks = 254, na.rm = T, dig.lab = -1)), '\\(|\\,|\\]'), function(m) m[[2]])))
        })
      } else {
        as.numeric(unlist(lapply(strsplit(levels(cut(i[i>=scalecut], breaks = 254, na.rm = T, dig.lab = -1)), '\\(|\\,|\\]'), function(m) m[[2]])))
      }
    })
    
    breakslist <- vector(mode = 'list', length = length(names(lowerlist)))
    names(breakslist) <- names(lowerlist)
    for (i in names(lowerlist)) {
      if (class(lowerlist[[i]])=='list') {
        breakslist[[i]] <- vector(mode = 'list', length = length(names(lowerlist[[i]])))
        names(breakslist[[i]]) <- names(lowerlist[[i]])
        for (j in names(lowerlist[[i]])) {
          breakslist[[i]][[j]] <- c(lowerlist[[i]][[j]], higherlist[[i]][[j]])
        }
      } else {
        breakslist[[i]] <- c(lowerlist[[i]], higherlist[[i]])
      }
    }
    # lower <- as.numeric(unlist(lapply(strsplit(levels(cut(mat[mat<scalecut], breaks = 2, na.rm = T, dig.lab = -1)), '\\(|\\,|\\]'), function(i) i[[2]])))
    # higher <- as.numeric(unlist(lapply(strsplit(levels(cut(mat[mat>=scalecut], breaks = 254, na.rm = T, dig.lab = -1)), '\\(|\\,|\\]'), function(i) i[[2]])))
    # breaks <- c(lower, higher)
    
    
    if (!is.null(mat2)) {
      tmp <- anncolslist[["consTOM"]]
      anncolslist[["consTOM"]] <- vector(mode = 'list', length = length(setLabels))
      names(anncolslist[['consTOM']]) <- setLabels
      for (i in setLabels) {
        anncolslist[['consTOM']][[i]] <- tmp
      }
    }
    for (i in names(anncolslist)) {
      if (length(anncolslist[[i]])==2) {
        for (j in names(anncolslist[[i]])) {
          anncolslist[[i]][[j]]$ModuleMembership <- anncolslist[[i]][[j]]$ModuleMembership[coltoplotlist[[i]][[j]]]
        }
      } else {
        anncolslist[[i]]$ModuleMembership <- anncolslist[[i]]$ModuleMembership[coltoplotlist[[i]]]
      }
    }
    
    # anncols$ModuleMembership <- anncols$ModuleMembership[coltoplot]
    
    for (i in names(matlist)) {
      if (class(matlist[[i]])=='list') {
        for (j in names(matlist[[i]])) {
          png(file.path(outpath, 'PLOTS', paste0('heatmap_colorder_', i, '_', j, '.png')), width = (210-32.5*2)*2, height = (210-32.5*2)*2, units = "mm", res = 300, pointsize = 10)
          pheatmap(matlist[[i]][[j]]$data[rowordlist[[i]], # [order(traits$LeafClassNumeric),]
                       anndtlist[[i]][[j]][col%in%coltoplotlist[[i]][[j]],rn]],
                   na_col = 'black',
                   # scale = "column",
                   # clustering_distance_rows = rowdist,
                   # clustering_distance_cols = coldist,
                   cluster_rows = F,
                   cluster_cols = F,
                   clustering_method = agglomeration,
                   show_colnames = F,
                   breaks = breakslist[[i]][[j]],
                   color = colorRampPalette(c('black','yellow'))(255),
                   annotation_row = rowannlist[[i]][[j]][,orderlists$Trait[[i]]$traitorder[orderlists$Trait[[i]]$traitorder%in%plottedTraits],drop=F],
                   annotation_col = colannlist[[i]][[j]][colannlist[[i]][[j]]$ModuleMembership%in%coltoplotlist[[i]][[j]],,drop=F],
                   annotation_colors = anncolslist[[i]][[j]])
          dev.off()
        }
      } else {
        png(file.path(outpath, 'PLOTS', paste0('heatmap_colorder_', i, '.png')), width = (210-32.5*2)*2, height = (210-32.5*2)*2, units = "mm", res = 300, pointsize = 10)
        pheatmap(matlist[[i]][rowordlist[[i]], # [order(traits$LeafClassNumeric),]
                                        anndtlist[[i]][col%in%coltoplotlist[[i]],rn]],
                 na_col = 'black',
                 # scale = "column",
                 # clustering_distance_rows = rowdist,
                 # clustering_distance_cols = coldist,
                 cluster_rows = F,
                 cluster_cols = F,
                 clustering_method = agglomeration,
                 show_colnames = F,
                 breaks = breakslist[[i]],
                 color = colorRampPalette(c('black','yellow'))(255),
                 annotation_row = rowannlist[[i]][,orderlists$Trait[[i]]$traitorder[orderlists$Trait[[i]]$traitorder%in%plottedTraits],drop=F],
                 annotation_col = colannlist[[i]][colannlist[[i]]$ModuleMembership%in%coltoplotlist[[i]],,drop=F],
                 annotation_colors = anncolslist[[i]])
        dev.off()
      }
    }
    # png(file.path(outpath, 'PLOTS', 'heatmap_colorder.png'), width = (210-32.5*2)*2, height = (210-32.5*2)*2, units = "mm", res = 300, pointsize = 10)
    # pheatmap(mat[roword, # [order(traits$LeafClassNumeric),]
    #              anndt[col%in%coltoplot,rn]],
    #          na_col = 'black',
    #          clustering_distance_rows = rowdist,
    #          clustering_distance_cols = coldist,
    #          cluster_rows = F,
    #          cluster_cols = F,
    #          clustering_method = agglomeration,
    #          show_colnames = F,
    #          breaks = breaks,
    #          color = colorRampPalette(c('black','yellow'))(255),
    #          annotation_row = rowann[,traitorder[traitorder%in%plottedTraits],drop=F],
    #          annotation_col = colann[colann$ModuleMembership%in%coltoplot,,drop=F],
    #          annotation_colors = anncols)
    # dev.off()
    # 
    # development = as.data.frame(traits$LeafClassNumeric)
    # names(development) = "LeafClassNumeric"
    # 
    # geneTraitSignificance = as.data.frame(cor(mat, development, use = "p"))
    # GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
    # # 
    # names(geneTraitSignificance) = paste("GS.", names(development), sep="")
    # names(GSPvalue) = paste("p.GS.", names(development), sep="")
    # # 
    if (!is.null(mat2)) {
      anndtlist$consTOM <- anndtlist$consTOM[[1]]
    }
  }
  
  
  # Create the starting data frame
  geneInfo0list <- vector(mode = 'list', length = length(matlist))
  names(geneInfo0list) <- names(matlist)
  for (i in names(matlist)) {
    if (class(matlist[[i]])=='list') {
      # geneInfo0list[[i]] <- vector(mode = 'list', length = length(matlist[[i]]))
      # names(geneInfo0list[[i]]) <- names(matlist[[i]])
      # for (j in names(matlist[[i]])) {
      geneInfo0list[[i]] <- data.table(id = colnames(matlist[[i]][[1]]$data),
                                       moduleColor = moduleColors[[i]])
      setkey(geneInfo0list[[i]], id)
      # }
    } else {
      geneInfo0list[[i]] <- data.table(id = colnames(matlist[[i]]),
                                            moduleColor = moduleColors[[i]])
      setkey(geneInfo0list[[i]], id)
    }
  }
  # geneInfo0 = data.table(id = colnames(mat),
  #                        moduleColor = moduleColors)
  
  
  # geneTraitSignificance,
  # GSPvalue)
  # setkey(geneInfo0, id)

  # HERE HERE HERE HERE HERE
  # for (i in names(anndtlist)) {
  #   if (all(class(anndtlist[[i]])=='list')) {
  #     for (j in names(anndtlist[[i]])) {
  #       anndtlist[[i]][[j]][,col:=factor(col, levels = sort(unique(col)), ordered = T, labels = sort(unique(col)))]
  #       setkey(anndtlist[[i]][[j]], col, ord)
  #     }
  #   } else {
  #     anndtlist[[i]][,col:=factor(col, levels = sort(unique(col)), ordered = T, labels = sort(unique(col)))]
  #     setkey(anndtlist[[i]], col, ord)
  #   }
  # }
  
  
  lapply(anndtlist, function(i) i[,col:=factor(col, levels = sort(unique(col)), ordered = T, labels = sort(unique(col)))])
  lapply(anndtlist, setkey, col, ord)
  
  # anndt[,col:=factor(col, levels = colord, ordered = T, labels = colord)]
  # setkey(anndt, col, ord)
  
  for (i in names(anndtlist)) {
    if (all(class(anndtlist[[i]])=='list')) {
      for (j in names(anndtlist[[i]])) {
        geneInfo0list[[i]][[j]] <- geneInfo0list[[i]][[j]][J(anndtlist[[i]][[j]][,rn]),]
      }
    } else {
      geneInfo0list[[i]] <- geneInfo0list[[i]][J(anndtlist[[i]][,rn]),]
    }
  }
  # geneInfo0 <- geneInfo0[.(anndt[,rn]),]
  
  lapply(names(geneInfo0list), function(i) fwrite(geneInfo0list[[i]], file = file.path(outpath, 'CSVs', paste0('wgcna_', i, '.csv')), sep = ',', col.names = T, row.names = F))
  # fwrite(geneInfo0, file = file.path(outpath, 'CSVs', 'wgcna.csv'), sep = ',', col.names = T, row.names = F)
  
  lapply(names(geneInfo0list), function(i) {
    dir.create(file.path(outpath, 'CSVs', 'Modules', i), recursive = T, showWarnings = F)
    geneInfo0list[[i]][,fwrite(.SD, file.path(outpath, 'CSVs', 'Modules', i, paste0(moduleColor, '.csv')), sep = ',', col.names = F, row.names = F),by=list(moduleColor)]
  })
  # geneInfo0[,fwrite(.SD, file.path(outpath, 'CSVs', paste0('module_', moduleColor, '.csv')), sep = ',', col.names = F, row.names = F),by=list(moduleColor)]
  
  universe <- readRDS(file.path(outpath, 'RDS', 'universe.rds'))
  fwrite(data.table(universe), file = file.path(outpath, 'CSVs', 'go_background.csv'), sep = ',', col.names = F, row.names = F)
  
  # cl <- makePSOCKcluster(3)
  if (any(runnew[c('identifier', 'softPower', 'MEDissThres', 'agglomeration', 'outpath', 'scale', 'quantcut', 'minModuleSize', 'cutoff')])) {
    cat('Annotation enrichment!\n')
    # doGO(universe = universe,
    #      geneInfo0 = geneInfo0,
    #      nodesize = nodesize,
    #      golist = golist,
    #      cutoff = cutoff,
    #      outpath = outpath)
    lapply(names(geneInfo0list), function(i) {
      cat(i, '\n')
      gores <- quickGO(grouptable = parseGI0(geneInfo0list[[i]], 'moduleColor', universe),
                       gsubid = gsubid,
                       subtypename='moduleColor',
                       func_ann = func_ann,
                       identifier = identifier,
                       alternative="two.sided",
                       method.subtype = "none",
                       method.all = 'BH',
                       digits = 3,
                       cutoff = cutoff,
                       multimembersep = NULL,
                       subtypestoexclude = 'background',
                       go_rollup = func_ann$go_offspring,
                       ncores = 20)
      gores <- rbindlist(lapply(names(gores), function(j) gores[[j]][,'GeneSet':=j]))
      fwrite(gores, file.path(outpath, 'CSVs', 'GO', paste0('goall_', i, '.csv')), sep = ',', row.names = F, col.names = T)  
    })
    # gores <- quickGO(grouptable = parseGI0(geneInfo0, 'moduleColor', universe),
    #                  gsubid = '_p[STY][0-9]*$',
    #                  subtypename='moduleColor',
    #                  func_ann = func_ann,
    #                  identifier = identifier,
    #                  alternative="two.sided",
    #                  method.subtype = "none",
    #                  method.all = 'BH',
    #                  digits = 3,
    #                  cutoff = cutoff,
    #                  multimembersep = NULL,
    #                  subtypestoexclude = 'background',
    #                  go_rollup = func_ann$go_offspring,
    #                  ncores = 20)
    # gores <- rbindlist(lapply(names(gores), function(i) gores[[i]][,'GeneSet':=i]))
    # fwrite(gores, file.path(outpath, 'CSVs', 'GO', 'goall.csv'), sep = ',', row.names = F, col.names = T)

    # Cytoscape export
    dir.create(file.path(outpath, 'CYT'), recursive = T, showWarnings = F)
    setkeyv(func_ann$annotations, identifier)
    cat('Read TOM!\n')
    if (!exists('TOMlist')) {
      TOMlist <- readRDS(file.path(outpath, 'RDS', 'TOMlist.rds'))
    }
    if (automatedtomcut) {
      tomcut <- lapply(TOMlist, findTOMCUT)
      tomcut <- findTOMCUT(TOM)
    } else {
      tomcut <- tomcut
    }
    # for (i in ModuleMembership) {
    require(pbmcapply)
    cat('Write modules to disk!\n')
    lapply(names(ModuleMembershiplist), function(j) {
      pbmclapply(ModuleMembershiplist[[j]], function(i) {
        # cat(i, '\n')
        dir.create(file.path(outpath, 'CYT', j, i), recursive = T, showWarnings = F)
        # Select module probes
        if (class(matlist[[j]])=='list') {
          probes = colnames(matlist[[j]][[1]]$data)
        } else {
          probes = colnames(matlist[[j]])
        }
        # probes = colnames(mat)
        inModule = is.finite(match(moduleColors[[j]], i))
        modProbes = probes[inModule]
        modGenes = func_ann$annotations[gsub('_p[STY][0-9]*$','', modProbes),SYMBOL,mult='first']
        # Select the corresponding Topological Overlap
        modTOM = TOMlist[[j]][inModule, inModule]
        dimnames(modTOM) = list(modProbes, modProbes)
        # Export the network into edge and node list files Cytoscape can read
        invisible(exportNetworkToCytoscape(modTOM,
                                       edgeFile = file.path(outpath, 'CYT', j, i, paste("CytoscapeInput-edges_", j, '_', i, ".txt", sep="")),
                                       nodeFile = file.path(outpath, 'CYT', j, i, paste("CytoscapeInput-nodes_", j, '_', i, ".txt", sep="")),
                                       weighted = TRUE,
                                       threshold = tomcut,
                                       nodeNames = modProbes,
                                       altNodeNames = modGenes,
                                       nodeAttr = moduleColors[[j]][inModule],
                                       includeColNames = T))
      }, mc.cores = 10)
    })
    # pbmclapply(ModuleMembership, function(i) {
    #   # cat(i, '\n')
    #   dir.create(file.path(outpath, 'CYT', i), recursive = T, showWarnings = F)
    #   # Select module probes
    #   probes = colnames(mat)
    #   inModule = is.finite(match(moduleColors, i))
    #   modProbes = probes[inModule]
    #   modGenes = func_ann$annotations[gsub('_p[STY][0-9]*$','', modProbes),SYMBOL,mult='first']
    #   # Select the corresponding Topological Overlap
    #   modTOM = TOM[inModule, inModule]
    #   dimnames(modTOM) = list(modProbes, modProbes)
    #   # Export the network into edge and node list files Cytoscape can read
    #   cyt = exportNetworkToCytoscape(modTOM,
    #                                  edgeFile = file.path(outpath, 'CYT', i, paste("CytoscapeInput-edges_", i, ".txt", sep="")),
    #                                  nodeFile = file.path(outpath, 'CYT', i, paste("CytoscapeInput-nodes_", i, ".txt", sep="")),
    #                                  weighted = TRUE,
    #                                  threshold = tomcut,
    #                                  nodeNames = modProbes,
    #                                  altNodeNames = modGenes,
    #                                  nodeAttr = moduleColors[inModule],
    #                                  includeColNames = T)
    # }, mc.cores = 10)
    # }
    if (writeentirenetwork) {
      lapply(names(TOMlist), function(j) {
        dir.create(file.path(outpath, 'CYT', j, 'ALL'), recursive = T, showWarnings = F)
        if (class(matlist[[j]])=='list') {
          probes = colnames(matlist[[j]][[1]]$data)
        } else {
          probes = colnames(matlist[[j]])
        }
        # probes = colnames(mat)
        genes = func_ann$annotations[probes,SYMBOL,mult='first']
        cat('Write network to disk!\n')
        # Export the network into edge and node list files Cytoscape can read
        cyt = exportNetworkToCytoscape(TOMlist[[j]],
                                       edgeFile = file.path(outpath, 'CYT', j, 'ALL', paste("CytoscapeInput-edges_", j, '_ALL', ".txt", sep="")),
                                       nodeFile = file.path(outpath, 'CYT', j, 'ALL', paste("CytoscapeInput-nodes_", j, '_ALL', ".txt", sep="")),
                                       weighted = TRUE,
                                       threshold = tomcut,
                                       nodeNames = probes,
                                       altNodeNames = genes,
                                       nodeAttr = moduleColors[[j]],
                                       includeColNames = T)  
      })
      # dir.create(file.path(outpath, 'CYT', 'ALL'), recursive = T, showWarnings = F)
      # probes = colnames(mat)
      # genes = func_ann$annotations[probes,SYMBOL,mult='first']
      # cat('Write network to disk!\n')
      # # Export the network into edge and node list files Cytoscape can read
      # cyt = exportNetworkToCytoscape(TOM,
      #                                edgeFile = file.path(outpath, 'CYT', 'ALL', paste("CytoscapeInput-edges_", 'ALL', ".txt", sep="")),
      #                                nodeFile = file.path(outpath, 'CYT', 'ALL', paste("CytoscapeInput-nodes_", 'ALL', ".txt", sep="")),
      #                                weighted = TRUE,
      #                                threshold = tomcut,
      #                                nodeNames = probes,
      #                                altNodeNames = genes,
      #                                nodeAttr = moduleColors,
      #                                includeColNames = T)
    }
    # cat('Combine network and correlation!\n')
    # ls()[!grepl('outpath', ls())]
    # source('/media/msdata5/users_files/Martin Frejno/postdoc/projects/arabidopsis/src/parseTOMsites.R')
    # parseTOMsites(outpath)
  }
  # stopCluster(cl = cl)
}
# 
# mat <- readRDS('preproc/leaf.rds')
# # golist <- readRDS('preproc/golist.rds')
# traits <- readRDS('preproc/traits.rds')
# func_ann <- readRDS('preproc/ath_ann.rds')
# 
# doWGCNA(mat = mat,
#         traits = traits,
#         func_ann = func_ann,
#         identifier = 'TAIR',
#         softPower = 10,
#         MEDissThres = 0.1,
#         agglomeration = 'ward.D2',
#         colord = 'sorted',
#         selectedTraits = c('LeafClassNumeric', '0', '1', '2', '3'),
#         plottedTraits = c('LeafClassNumeric'),
#         outpath = 'res/20171011_wgcna',
#         scale = T,
#         scalecut = 0,
#         sampleord = rownames(traits),
#         quantcut = 0.75,
#         minModuleSize = 30,
#         nodesize = 10,
#         cutoff = 0.05)
