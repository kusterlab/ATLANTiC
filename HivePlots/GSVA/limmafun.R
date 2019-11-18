limmafun <- function(mat, design, p.adjust = 'BH', cutoff = 0.05, allpairs = F) {
  require(limma)
  require(Biobase)
  require(data.table)
  data <- mat[,names(design)]
  dataexp <- new("ExpressionSet", exprs=data)
  
  designmatrix <- model.matrix(~ -1+factor(design))
  rownames(designmatrix) <- names(design)
  colnames(designmatrix) <- levels(factor(design))
  
  if (allpairs) {
    # all pairwise comparisons
    z <- t(combn(colnames(designmatrix),m=2))
    con <- paste(z[,1], z[,2], sep="-")
    # Generate contrast.matrix
    contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(con),levels=list(designmatrix))))
  } else {
    confun <- function(i,vec) {
      paste0(i, '-(', paste0(vec[-match(i,vec)],collapse='+'), ')/', length(vec)-1)
    }
    con <- sapply(colnames(designmatrix), function(i) confun(i, colnames(designmatrix)))
    # Generate contrast.matrix
    contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(con),levels=list(designmatrix))))
  }
  fit <- lmFit(dataexp, designmatrix)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  efit <- eBayes(fit2)
  tab <- data.table(topTable(efit, adjust.method=p.adjust, number="all", p.value = 1), keep.rownames=T)
  setkey(tab, rn)
  tab2 <- melt(tab[,c('rn', unique(design)),with=F], id.vars = 'rn', value.name = 'log2-fold-change')
  tab <- tab[adj.P.Val < cutoff,]
  res <- decideTests(object = efit, method = "separate", adjust.method = "BH", p.value = cutoff)
  sigs <- apply(res@.Data, 2, function(i) names(which(i==1)))
  # efit$p.value <- matrix(p.adjust(efit$p.value, method = p.adjust), nrow = nrow(efit$p.value), byrow = F, dimnames = dimnames(efit$p.value))
  
  # efit$p.value <- apply(efit$p.value, 2, function(i) p.adjust(i, method = p.adjust))
  efit$p.value[abs(res)!=1] <- NA
  res_dt <- data.table(melt(res))
  setnames(res_dt, c('Variable', 'Group', 'Direction'))
  res_dt <- res_dt[Direction!=0]
  res_dt[,Variable:=as.character(Variable)]
  res_dt[,Group:=as.character(Group)]
  res_dt[,Direction:=as.character(c('down', '', 'up')[Direction+2])]
  setkey(res_dt, Variable, Group)
  pvals_dt <- data.table(melt(as.matrix(efit$p.value)))
  setnames(pvals_dt, c('Variable', 'Group', 'adj.P.Val'))
  pvals_dt <- pvals_dt[!is.na(adj.P.Val)]
  pvals_dt[,Variable:=as.character(Variable)]
  pvals_dt[,Group:=as.character(Group)]
  setkey(pvals_dt, Variable, Group)
  pvals_dt <- pvals_dt[res_dt][,list(Variable, Group, Direction, adj.P.Val)]
  setkey(pvals_dt, Variable, Group)
  setkey(tab2, rn, variable)
  pvals_dt[tab2,`log2-fold-change`:=`log2-fold-change`]
  # combs <- combn(names(sigs), m = 2)
  # ints <- apply(combs, 2, function(i) intersect(sigs[[i[1]]],sigs[[i[2]]]))
  # names(ints) <- apply(combs, 2, function(i) paste0(i, collapse = '-'))
  # tab[,plot(orange3_vs_not_orange3,-log10(adj.P.Val))]
  # tab[.(sigs$orange3_vs_not_orange3),points(orange3_vs_not_orange3,-log10(adj.P.Val), col=3)]
  # tab[.(ints$`orange3_vs_not_orange3-red3_vs_not_red3`),points(orange3_vs_not_orange3,-log10(adj.P.Val), col=2)]
  # boxplot(mat[ints$`green4_vs_not_green4-magenta4_vs_not_magenta4`,names(design)], col = design, las = 2)
  return(list(tab=tab,
              pvals=pvals_dt,
              res=res,
              sigs=sigs))
}

parselim_o <- function(lim) {
  require(data.table)
  setkey(lim$tab, rn)
  res <- rbindlist(lapply(names(lim$sigs), function(i) lim$tab[.(lim$sigs[[i]]),list(source=rn,target=i,weight=adj.P.Val),nomatch = 0])) # ,logFC=eval(as.name(i))
  return(res)
}

parselim <- function(lim) {
  require(data.table)
  setkey(lim$pvals, Direction)
  res <- lim$pvals[.('up')]
  res[,Direction:=NULL]
  setnames(res, c('source', 'target', 'weight', 'FC'))
  return(res)
}

filterLimmaRes_old <- function(lim, direction) {
  direction <- match.arg(direction, choices = c('up', 'down', 'both'))
  if (direction=='both')
    direction <- c('up', 'down')
  setkey(lim$tab, rn)
  pos <- apply(lim$res@.Data, 2, function(i) names(which(i==1)))
  neg <- apply(lim$res@.Data, 2, function(i) names(which(i==-1)))
  sitespos <- rbindlist(lapply(names(pos), function(i) lim$tab[.(pos[[i]]),list(Variable=rn, Group=i, Direction = 'up', log2fc=eval(as.name(i)), adj.P.Val),nomatch = 0]))
  sitesneg <- rbindlist(lapply(names(neg), function(i) lim$tab[.(neg[[i]]),list(Variable=rn, Group=i, Direction = 'down', log2fc=eval(as.name(i)), adj.P.Val),nomatch = 0]))
  result <- rbindlist(list(sitespos, sitesneg))
  setkey(result, Direction)
  result <- result[.(direction),]
  return(result)
}

filterLimmaRes <- function(lim, direction) {
  direction <- match.arg(direction, choices = c('up', 'down', 'both'))
  if (direction=='both')
    direction <- c('up', 'down')
  setkey(lim$tab, rn)
  setkey(lim$pvals, Variable, Group, Direction)
  pos <- apply(lim$res@.Data, 2, function(i) names(which(i==1)))
  neg <- apply(lim$res@.Data, 2, function(i) names(which(i==-1)))
  sitespos <- rbindlist(lapply(names(pos), function(i) lim$tab[.(pos[[i]]),list(Variable=rn, Group=i, Direction = 'up', log2fc=eval(as.name(i))),nomatch = 0]))
  setkey(sitespos, Variable, Group, Direction)
  sitespos <- sitespos[lim$pvals,nomatch=0]
  sitesneg <- rbindlist(lapply(names(neg), function(i) lim$tab[.(neg[[i]]),list(Variable=rn, Group=i, Direction = 'down', log2fc=eval(as.name(i))),nomatch = 0]))
  setkey(sitesneg, Variable, Group, Direction)
  sitesneg <- sitesneg[lim$pvals,nomatch=0]
  result <- rbindlist(list(sitespos, sitesneg))
  setkey(result, Direction)
  result <- result[.(direction),]
  return(result)
}

impute <- function(x, imputation="rawfilemin", columnwise=F, downshift=1.8, width=0.3, mc.cores = 10, seed = 23021987, evmap_red = NULL, timerange_min = 5) {
  RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = "default")
  set.seed(seed)
  if (imputation=="rawfilemin") {
    imputed <- x
    rownames(imputed) <- paste0(rownames(imputed), "_", 1:nrow(imputed))
    for (i in 1:dim(imputed)[2]) {
      mini <- min(imputed[,i], na.rm=T)
      imputed[is.na(imputed[,i]),i] <- mini
    }
    rownames(imputed) <- gsub("_[0-9]*", "", rownames(imputed))
  } else {
    if (imputation=="perseus") {
      imputed <- x
      rownames(imputed) <- paste0(rownames(imputed), "_", 1:nrow(imputed))
      if (columnwise==F) {
        meanx <- mean(x, na.rm=T)
        sdx <- sd(x, na.rm=T)
        meanimpute <- meanx-sdx*downshift
        sdimpute <- sdx*width
        imputed[is.na(imputed)] <- rnorm(n=length(which(is.na(x))), mean=meanimpute, sd=sdimpute)
      } else {
        meanx <- apply(x, 2, function(dat) mean(dat, na.rm=T))
        sdx <- apply(x, 2, function(dat) sd(dat, na.rm=T))
        meanimpute <- meanx-sdx*downshift
        sdimpute <- sdx*width
        for (i in 1:dim(imputed)[2]) {
          imputed[is.na(imputed[,i]),i] <- rnorm(n=length(which(is.na(x[,i]))), mean=meanimpute[i], sd=sdimpute[i])
        }
      }
      rownames(imputed) <- gsub("_[0-9]*", "", rownames(imputed))
    } else {
      if (imputation=="proteinmin") {
        colidx <- sapply(groupnames, function(dat) grep(paste("\\b[", paste(strsplit(dat, split="")[[1]], collapse="]["), "]\\b", sep=""), colnames(x)))
        imputed <- matrix(ncol=ncol(x), nrow=nrow(x))
        for (i in 1:dim(x)[1]) {
          # Protein-wise minimum imputation
          impute <- apply(colidx, 2, function(dat) x[i,dat])
          kinmin <- min(impute, na.rm=T)
          impute[which(is.na(impute))] <- kinmin
          longimpute <- unlist(apply(impute, 1, function(dat) list(dat)))
          imputed[i,] <- longimpute}
      } else {
        if (imputation=="smallerthanx") {
          # set.seed(23021987)
          # smallerthanx imputation as discussed
          require(data.table)
          require(reshape2)
          print("Imputation not possible column-wise!")
          imputed <- x
          rownames(imputed) <- paste0(rownames(imputed), "_", 1:nrow(imputed))
          longx <- data.table(melt(imputed))
          setkey(longx, Var1)
          amounts <- longx[is.na(value),length(unique(Var2)), by=Var1]
          amounts[,Var1:=as.character(Var1)]
          setkey(amounts, Var1)
          setkey(longx, Var1)
          mins <- longx[amounts[,Var1],min(value, na.rm=T),by=Var1]
          mins[,Var1:=as.character(Var1)]
          mins <- mins[is.finite(V1),]
          onlyNAs <- setdiff(amounts$Var1, mins$Var1)
          amounts <- amounts[Var1%in%mins[,Var1],]
          setkey(mins, Var1)
          setkey(longx, value)
          longxnona <- longx[!is.na(value),]
          setkey(longx, Var1)
          
          # new version - lower RAM & data.table output
          mc.reset.stream()
          repl <- data.table(Var1=as.character(longx[is.na(value)][amounts[,Var1],Var1]), Var2=as.character(longx[is.na(value)][amounts[,Var1],Var2]), rbindlist(mclapply(as.character(amounts[,Var1]), function(i) longxnona[seq_len(longxnona[J(mins[i,V1]), which=TRUE, mult="first"])][,list(value=sample(value, size = amounts[i,V1], replace = T))], mc.cores=mc.cores, mc.set.seed = F)))
          setkey(repl, Var1)
          repl[as.character(longx[which.min(value),Var1]),value:=longx[which.min(value),value]]
          longx <- rbind(longxnona, repl)
          
          #         
          # previous version - lower RAM, but still list
          #           repl <- mclapply(as.character(amounts[,Var1]), function(i) longxnona[seq_len(longxnona[J(mins[i,V1]), which=TRUE, mult="first"])][,sample(value, size = amounts[i,V1], replace = T)], mc.cores=mc.cores)
          #           names(repl) <- amounts[,Var1]
          #           if (!is.null(repl[[as.character(longxnona[which.min(value),Var1])]])) {
          #             repl[[as.character(longxnona[which.min(value),Var1])]] <- rep(longxnona[which.min(value),value], times=amounts[as.character(longxnona[which.min(value),Var1]),V1])
          #           }
          #           urepl <- unlist(repl)
          #           longx[is.na(value),value:=urepl]
          
          tmp <- acast(data = longx, formula = Var1~Var2, value.var = "value")
          tmp <- rbind(tmp, imputed[onlyNAs,])
          imputed <- tmp[rownames(imputed),colnames(imputed)]
          rownames(imputed) <- gsub("_[0-9]*", "", rownames(imputed))
        } else {
          if (imputation=="mix") {
            # step-wise knn and smallerthanx imputation
            require(data.table)
            require(reshape2)
            print("Imputation not possible column-wise!")
            imputed <- eval(x)
            rownames(imputed) <- paste0(rownames(imputed), "_", 1:nrow(imputed))
            nullidx <- imputed[imputed==0]
            imputed_iknn <- impute.knn(data = imputed, rng.seed = seed, k = 10)$data
            remidy <- sapply(which(apply(imputed, 1, function(i) any(is.na(i)))), function(i) intersect(which(alldup(imputed_iknn[i,])),which(is.na(imputed[i,]))))
            remidy <- remidy[which(lapply(remidy, length)>0)]
            for (i in seq_along(remidy)) {
              imputed_iknn[names(remidy[i]),remidy[[i]]] <- NA
            }
            remidx <- sapply(1:ncol(imputed_iknn), function(i) intersect(which(alldup(imputed_iknn[,i])),which(is.na(imputed[,i]))))
            names(remidx) <- colnames(imputed_iknn)
            for (i in seq_along(remidx)) {
              imputed_iknn[remidx[[i]],i] <- NA
            }
            imputed_iknn[imputed_iknn==0] <- NA
            imputed_iknn[nullidx] <- 0
            rownames(imputed_iknn) <- gsub("_[0-9]*", "", rownames(imputed_iknn))
            cat(paste0(as.character(round(1-length(which(is.na(imputed_iknn)))/length(which(is.na(imputed))), digits = 3)*100), "%"), "of missing values imputed with KNN!")
            # set.seed(23021987)
            imputed_imix <- log2(impute2(2^imputed_iknn, imputation = "smallerthanx", columnwise = F, mc.cores = mc.cores))
            rownames(imputed_imix) <- gsub("_[0-9]*", "", rownames(imputed_imix))
            imputed <- list(iknn=imputed_iknn, imix=imputed_imix)
          } else {
            if (imputation=="truncnorm") {
              require(truncnorm)
              require(mclust)
              imputed <- x
              rownames(imputed) <- paste0(rownames(imputed), "_", 1:nrow(imputed))
              if (columnwise) {
                # medians <- apply(imputed, 1, function(i) median(i, na.rm=T))
                mins <- apply(imputed, 1, function(i) min(i, na.rm=T))
                distlist <- apply(imputed, 2, function(i) mins[which(is.na(i))])
                mlist <- lapply(distlist, function(i) densityMclust(i, G=1))
                # hist(rtruncnorm(n=length(test[[1]]), a = -Inf, b = 15, mean = m1$parameters$mean, sd = sqrt(m1$parameters$variance$sigmasq)), nclass=100)
                mc.reset.stream()
                for (i in 1:ncol(imputed)) {
                  imputed[is.na(imputed[,i]),i] <- unlist(mclapply(mins[is.na(imputed[,i])], function(j) rtruncnorm(n=1, a = min(distlist[[i]], na.rm=T), b = j, mean = mlist[[i]]$parameters$mean, sd = sqrt(mlist[[i]]$parameters$variance$sigmasq)), mc.cores = mc.cores, mc.set.seed = F))
                }
              } else {
                narows <- which(rowSums(is.na(imputed))>0)
                mins <- apply(imputed, 1, function(i) min(i, na.rm=T))
                dist <- unlist(apply(imputed, 2, function(i) mins[which(is.na(i))]))
                model <- densityMclust(dist, G=1)
                for (i in names(narows)) {
                  imputed[i,is.na(imputed[i,])] <- rtruncnorm(n=length(imputed[i,is.na(imputed[i,])]), a = min(dist, na.rm=T), b = mins[i], mean = model$parameters$mean, sd = sqrt(model$parameters$variance$sigmasq))
                }
                imputed[which.min(mins),is.na(imputed[which.min(mins),])] <- rep(min(mins), times=length(imputed[which.min(mins),is.na(imputed[which.min(mins),])]))
              }
              rownames(imputed) <- gsub("_[0-9]*", "", rownames(imputed))
            } else {
              if (imputation=="pepint") {
                stopifnot(!is.null(evmap_red))
                print("Preparing...")
                imputed <- x
                nagenes <- names(which(rowSums(is.na(imputed))>0)) # Genes with at least one NA
                setkey(evmap_red, EXPID, GENEID, iBAQ_gn_mc)
                protint <- evmap_red[colnames(imputed),sum(iBAQ_gn_mc, na.rm=T),by=list(EXPID, GENEID)] # Prot Quant based on evmap_red
                invisible(protint[V1==0,V1:=NA])
                invisible(protint[,V1:=signif(log2(V1), 5)])
                protint <- protint[!is.na(V1)]
                setkey(protint, GENEID, V1) # sort by gene names AND THEN expression
                protint <- protint[nagenes,mult="first",nomatch=0] # reduce PROTEIN expression table to minimum expression ACROSS ALL SAMPLES
                setkey(protint, GENEID)
                mins <- evmap_red[.(protint[,EXPID],protint[,GENEID]),list(EXPID, GENEID, FRACID, RT, CRT, iBAQ_gn_mc)] # find minimum FEATURE expression making up minimum protein expression
                # plot(mins[,signif(log2(sum(iBAQ_gn_mc)),5),by=list(GENEID)][order(GENEID),V1],protint[,V1], pch='.')
                # lengthmins <- mins[,length(unique(RT)),by=list(FRACID, GENEID)] # count number of features per gene and fraction => multiple features (and therefore intensities) for the same gene in the same fraction (different peptides and different charge states of the same peptide)!
                # setkey(lengthmins, FRACID, GENEID)
                # invisible(mins[,id:=paste0(FRACID, '_|_', GENEID)])
                # setkey(mins, id, iBAQ_gn_mc)
                # mins <- mins[unique(mins$id),mult="first"]
                setkey(mins, FRACID, GENEID, iBAQ_gn_mc)
                setkey(evmap_red, FRACID, iBAQ_gn_mc)
                #                 sampdt <- rbindlist(lapply(mins[,unique(FRACID)], function(i) {
                #                   rbindlist(mclapply(mins[.(i),unique(GENEID)], function(j) {
                #                     data.table(FRACID=i, GENEID=j, EXPID=rep(names(which(is.na(imputed[j,]))), times=lengthmins[.(i,j),V1]), int=evmap_red[.(i),list(iBAQ_gn_mc)][iBAQ_gn_mc<=mins[.(i,j),iBAQ_gn_mc],sample(x = iBAQ_gn_mc, size = length(which(is.na(imputed[j,])))*lengthmins[.(i,j),V1], replace = T)])
                #                   }, mc.cores=mc.cores))
                #                 }))
                #                 sampdt <- rbindlist(lapply(mins[,unique(FRACID)], function(i) {
                #                   rbindlist(mclapply(mins[.(i),unique(GENEID)], function(j) {
                #                     rbindlist(lapply(mins[.(i,j),unique(iBAQ_gn_mc)], function(k) {
                #                       data.table(FRACID=i, GENEID=j, EXPID=names(which(is.na(imputed[j,]))), int=evmap_red[.(i),list(iBAQ_gn_mc)][!is.na(iBAQ_gn_mc) & iBAQ_gn_mc<=mins[.(i,j,k),iBAQ_gn_mc],sample(x = iBAQ_gn_mc, size = length(which(is.na(imputed[j,]))), replace = T)])
                #                     }))
                #                   }, mc.cores=mc.cores))
                #                 }))
                sampdt <- rbindlist(lapply(mins[,unique(FRACID)], function(i) {
                  print(paste0("Processing Fraction ", i, "..."))
                  mc.reset.stream()
                  rbindlist(mclapply(mins[.(i),unique(GENEID)], function(j) {
                    rbindlist(lapply(mins[.(i,j),unique(iBAQ_gn_mc)], function(k) {
                      data.table(FRACID=i, GENEID=j, EXPID=names(which(is.na(imputed[j,]))), int=evmap_red[.(i),list(CRT,iBAQ_gn_mc)][!is.na(iBAQ_gn_mc) & iBAQ_gn_mc<=k & CRT >= mins[.(i,j,k),CRT]-timerange_min/2 & CRT <= mins[.(i,j,k),CRT]+timerange_min/2,sample(x = iBAQ_gn_mc, size = length(which(is.na(imputed[j,]))), replace = T)])
                    }))
                  }, mc.cores = mc.cores, mc.set.seed = F))
                }))
                sampdt_genesum <- sampdt[,signif(log2(sum(int, na.rm=T)), 5), by=list(EXPID, GENEID)]
                setkey(sampdt_genesum, GENEID, EXPID)
                longx <- data.table(melt(imputed))
                setkey(longx, Var1, Var2)
                invisible(longx[sampdt_genesum,value:=V1])
                tmp <- acast(data = longx, formula = Var1~Var2, value.var = "value")
                imputed <- tmp[rownames(imputed),colnames(imputed)]
                #                 test <- rbindlist(mclapply(colnames(imputed)[1:2], function(i) {
                #                   # print(i)
                #                   genes <- rownames(imputed)[is.na(imputed[,i])]
                #                   tmplist <- list()
                #                   for (j in genes) {
                #                     # print(j)
                #                     fracs <- mins[j,FRACID]
                #                     for (k in fracs) {
                #                       tmplist[[j]][[k]] <- sample(evmap_red[.(k)][iBAQ_gn_mc<=mins[.(j,k),iBAQ_gn_mc],iBAQ_gn_mc], size = 1, replace = F)
                #                       tmplist[[j]] <- sum(tmplist[[i]][[j]], na.rm=T)
                #                     }
                #                   }
                #                   tmp <- data.table(sample=i, `Gene names`=names(tmplist), iBAQ_gn_mc=unlist(tmplist))
                #                 }, mc.cores = mc.cores))
              }
            }
          }
        }
      }
    }
  }
  if (is.list(imputed)) {
    for (i in seq_along(imputed)) {
      colnames(imputed[[i]]) <- colnames(x)
      rownames(imputed[[i]]) <- rownames(x)
    }
  } else {
    colnames(imputed) <- colnames(x)
    rownames(imputed) <- rownames(x)
  }
  print("Done")
  RNGkind(kind = "default", normal.kind = "default")
  return(imputed)
}
