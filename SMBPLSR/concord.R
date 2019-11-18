# norm <- function(x) x/c(sd(x))/sqrt(length(x)-1)
norm <- function(x) x/sqrt(c(crossprod(x)))

tnorm <- function(x, f) {
  v <- unique(f)
  newv <- rep(NA, length(x))
  for (i in v) {
    ii <- f == i
    xx <- norm(x[ii])
    newv[ii] <- xx
  }
  newv
}

concord <- function(X, Y, nf=2, dmod = 1) {
  
  nmat <- length(X)
  if (is.null(names(X)))
    names(X) <- paste0("X", 1:nmat)
  
  nr <- sapply(X, nrow)
  nc <- unique(sapply(X, ncol))
  if (length(nc) > 1)
    stop("Number of columns in X need to be the same.")
  
  i.sample <- rep(names(X), each = nc)
  i.feature <- rep(names(X), nr)
  
  Ynorm <- scale(t(Y), center = TRUE, scale = TRUE)
  Xnorm <- lapply(X, function(x) scale(t(x), center = TRUE, scale = TRUE))
  Xcat <- do.call("cbind", Xnorm)
  Ynorm.o <- Ynorm
  Xnorm.o <- Xnorm
  
  ys <- yloading <- gls <- bls <- loading <- c()
  xvar <- c()
  var <- c()

  for (f in 1:nf) {
    print(f)
    if (f == 1 || dmod != 1 )
      S <- t(Ynorm) %*% Xcat
    if (f == 1)
      S.o <- S
    
    decom <- propack.svd(S, neig = 1, opts = list(kmax = 20))
    
    xa <- Xcat %*% decom$v[, 1]
    yb <- Ynorm %*% decom$u[, 1]
    var <- c(var, decom$d[1]^2)
    ys <- cbind(ys, yb)
    gls <- cbind(gls, xa)
    yloading <- cbind(yloading, decom$u[, 1])
    loading <- cbind(loading, decom$v[, 1])
    
    xai <- lapply(names(X), function(x) {
      ii <- i.feature == x
      Xcat[, ii] %*% decom$v[ii, 1]
      })
    
    xai.var <- sapply(xai, crossprod)
    xvar <- cbind(xvar, xai.var)
    bls <- cbind(bls, unlist(xai))
    
    if (dmod == 1) {
      # deflation of S, crossprod matrix, rewrite do not use for loop
      # the classical concordance approach
      S <- S - tcrossprod(decom$u[, 1]) %*% S
      # or the same 
      # Ynorm <- Ynorm - Ynorm %*% tcrossprod(decom$u[, 1])
    } else if (dmod == 2) {
      # defaltion Y using its normed score
      Ynorm <- Ynorm - t(t(Ynorm) %*% tcrossprod(norm(yb)))
    } else if (dmod == 3) {
      # defaltion X and Y using normed score of Y
      Ynorm <- Ynorm - t(t(Ynorm) %*% tcrossprod(norm(yb)))
      Xcat <- Xcat - t(t(Xcat) %*% tcrossprod(norm(yb)))
    } else if (dmod == 4) {
      # deflaltion X using loading of X
      Xcat <- Xcat - Xcat %*% tcrossprod(decom$v[, 1])
    } else {
      stop("unknown dmod")
    }

  }
  
  rownames(loading) <- colnames(Xcat)
  
  list(i.sample = i.sample,
       i.feature = i.feature, 
       ys = ys, 
       ys.norm = apply(ys, 2, norm),
       yloading = yloading, 
       gls = gls,
       gls.norm = apply(gls, 2, norm),
       bls = bls,
       bls.norm = apply(bls, 2, tnorm, f=i.sample),
       loading = loading,
       xvar = xvar, 
       var = var, 
       S = S.o, 
       Ynorm = Ynorm.o, 
       Xnorm = Xnorm.o)
}



sconcord.mpbca <- function(X, Y, nf = 2, k, deflat = "globalScore") {
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
                 k = k, 
                 center = FALSE, 
                 scale = FALSE, 
                 option = "uniform",
                 method = deflat, 
                 svd.solver = "propack", 
                 moa = FALSE, 
                 verbose = FALSE)
    
    # block score of Xs
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
    Xnorm <- lapply(1:length(Xnorm), function(i) {
      Xv <- Xsv[, i, drop = FALSE]
      Xnorm[[i]] - Xv %*% t(Xv) %*% Xnorm[[i]]
    })
    names(Xnorm) <- Xnames
    
    # assign results
    res$loadingXs = cbind(res$loadingXs, unlist(mr1$pb))
    res$scoreXs = cbind(res$scoreXs, c(v1))
    res$scoreXall = cbind(res$scoreXall, v1 %*% (1/mr1$w))
    res$scoreXs1 = cbind(res$scoreXs, c(scale(v1, center = FALSE, scale = TRUE)))
    res$scoreXall1 = cbind(res$scoreXall, scale(res$scoreXall, center = FALSE, scale = TRUE))
    res$loadingY = cbind(res$loadingY, mr1$t)
    res$loadingY0 = cbind(res$loadingY0, sto/sqrt(c(crossprod(mr1$t))))
    res$scoreY = cbind(res$scoreY, tsparse)
    res$scoreY1 = cbind(res$scoreY1, scale(res$scoreY, center = FALSE, scale = TRUE))
    res$w <- cbind(res$w, mr1$w)
  } 
  res
  
}

checkCC <- function(x, tol = 1e-7) {
  orth <- function(x, t = tol) {
    all(x[upper.tri(x)] < t) & all(x[lower.tri(x)] < t)
  }
  ii <- x$var > tol
  l <- list(y.score = crossprod(x$ys[, ii]),
       y.score.norm = crossprod(x$ys.norm[, ii]), 
       x.score = crossprod(x$gls[, ii]),
       x.score.norm = crossprod(x$gls.norm[, ii]), 
       x.block.score = crossprod(x$bls[, ii]), 
       x.block.score.norm = crossprod(x$bls.norm[, ii]), 
       x.loading = crossprod(x$loading[, ii]), 
       y.loading = crossprod(x$yloading[, ii]), 
       x.y.score = crossprod(x$gls[, ii], x$ys[, ii]),
       x.y.score.cor = cor(x$gls[, ii], x$ys[, ii])
       )
  sapply(l, orth)
}
