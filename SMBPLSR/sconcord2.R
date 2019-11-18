# incomplete implementation of sconcord methods. 
# It's an extension of sconcord.mbpca function


sconcord2 <- function(X, Y, nf = 2, signif = 0.05, nperm = 30, deflat = "globalScore") {
  Ynorm <- scale(t(Y), center = TRUE, scale = TRUE)
  Xnorm <- lapply(X, function(x) scale(t(x), center = TRUE, scale = TRUE))
  Xnames <- names(Xnorm)
  
  res <- list(loadingXs = c(), 
              scoreXs = c(), 
              scoreXall = c(),
              scoreXs1 = c(),
              scoreXall1 = c(),
              loadingY = c(),
              scoreY = c(), 
              scoreY1 = c()
  )
  
  for (i in 1:nf) {
    cat(paste("calculating component", i, "...\n"))
    
    xl <- lapply(Xnorm, function(x) t(crossprod(Ynorm, x)))
    mr1 <- mbpca(xl, 
                 ncomp = 1, 
                 k = "all", 
                 center = FALSE, 
                 scale = FALSE, 
                 option = "uniform",
                 method = deflat, 
                 svd.solver = "propack", 
                 moa = FALSE, 
                 verbose = FALSE)
    
    rep <- replicate(nperm, {
      xlp <- lapply(Xnorm, function(x) t(crossprod(Ynorm[sample(1:nrow(Ynorm), replace = TRUE), ], x)))
      mr1p <- mbpca(xlp, ncomp = 1, k = "all", center = FALSE, scale = FALSE, option = "uniform",
                    method = deflat, svd.solver = "propack", moa = FALSE, verbose = FALSE)
      unlist(mr1p$pb)
    })
    
    spb <- unlist(res$pb)
    pm <- pmin(pnorm(q = abs(spb), mean = rowMeans(rep), sd = rowMads(rep), lower.tail = FALSE)*2, 1)
    
    
    
    res <- list(mr1, rep)
    
    ##################
    v1 <- mapply(function(x, y) c(x %*% y),
                 x = Xnorm, y = mr1$pb, SIMPLIFY = TRUE)

    fg <- function(r, n, cut = 0.005, plot = FALSE) {
      br <- (r[-1] - r[-length(r)])/r[-length(r)]
      pos <- max(which(br > cut))
      ns <- n[pos]
      if (plot) {
        layout(matrix(1:2, 1, 2))
        barplot(br)
        abline(h=cut)
        barplot(r)
        abline(h = r[ns])
      }
      ns
    }

    nss <- round(seq(1, ncol(Ynorm)-1, length.out = 100))
    corsum <- sapply(nss, function(ns) {
      v2 <- sapply(mr1$tb, function(x) {
        ist <- softK(x, ns)
        c(Ynorm %*% ist)
      })
      sum(diag(cor(v1, v2)))
    })
    sto <- softK(mr1$t, fg(corsum, nss))
    tsparse <-  Ynorm %*% sto
    mr1$corsum <- corsum

    Ysv <- norm(tsparse)
    Xsv <- apply(v1, 2, norm)
    Ynorm  <- Ynorm - Ysv %*% t(Ysv) %*% Ynorm
    Xnorm <- lapply(1:3, function(i) {
      Xv <- Xsv[, i, drop = FALSE]
      Xnorm[[i]] - Xv %*% t(Xv) %*% Xnorm[[i]]
    })
    names(Xnorm) <- Xnames
  #   
  #   res$loadingXs = cbind(res$loadingXs, unlist(mr1$pb))
  #   res$scoreXs = cbind(res$scoreXs, c(v1))
  #   res$scoreXall = cbind(res$scoreXall, v1 %*% (1/mr1$w))
  #   res$scoreXs1 = cbind(res$scoreXs, c(scale(v1, center = FALSE, scale = TRUE)))
  #   res$scoreXall1 = cbind(res$scoreXall, scale(res$scoreXall, center = FALSE, scale = TRUE))
  #   res$loadingY = cbind(res$loadingY, mr1$t)
  #   res$loadingY0 = cbind(res$loadingY0, sto/sqrt(c(crossprod(mr1$t))))
  #   res$scoreY = cbind(res$scoreY, tsparse)
  #   res$scoreY1 = cbind(res$scoreY1, scale(res$scoreY, center = FALSE, scale = TRUE))
  #   res$w <- cbind(res$w, mr1$w)
  }
  res
  # 
}
