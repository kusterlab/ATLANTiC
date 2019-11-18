library(RUnit)
#' The function returns
#' @title percentage of true values according to a rankable list
#' @param x numeric vector
#' @param f the logical vector
#' @param perc minimum percentage accepted
#' 
#' @return 
#' 
#' @example 
#' 
percTrue <- function(x, f, perc = 0.8) {
  if (any(x < 0))
    stop("only works for positive numbers")
  rk <- rank(-x)
  urk <- sort(unique(rk))
  acc <- sapply(urk, function(i) {
    ii <- rk <= i
    sum(f[ii])/length(f[ii])
  })
  ss <- which(acc > perc)
  if (length(ss) == 0)
    return(0)
  ceiling(urk[max(ss)])
}

test.percTrue <- function() {
  i <- sample(1:10)
  x <- rev(10*(1:10))[i]
  f <- rep(c(TRUE, FALSE), c(3, 7))[i]
  checkEqualsNumeric(percTrue(x, f, perc = 0.8), 3)
  checkEqualsNumeric(percTrue(x, f, perc = 0.7), 4)
  x <- 10*(1:10)
  f <- rep(c(TRUE, FALSE), c(3, 7))
  checkEquals(percTrue(x, f, perc = 0.7), 0)
}

#' 
#' @param x
#' @param npos
#' @param nneg
#' 
#' @return 
#' @example 

bisoft <- function(x, npos, nneg) {
  neg <- x <= 0
  
  # nneg <- min(nneg, sum(neg))
  # npos <- min(npos, sum(!neg))
  # 
  tv <- ifelse(sum(neg) <= nneg, 0, sort(x, decreasing = TRUE)[npos+1])
  bv <- ifelse(sum(!neg) <= npos, 0, sort(x)[nneg+1])
  
  x[neg] <- pmin(x[neg] - bv, 0)
  x[!neg] <- pmax(x[!neg] - tv, 0)
  c(x)
}

test.bisoft <- function() {
  x <- -10:10
  r <- bisoft(x, 4, 2)
  
  ii <- sample(1:21)
  xs <- x[ii]
  rs <- bisoft(xs, 4, 2)
  checkEqualsNumeric(rs, r[ii])
  
  rex <- bisoft(x, 10, 10)
  checkEqualsNumeric(x, rex)
  rex <- bisoft(x, 20, 20)
  checkEqualsNumeric(x, rex)
  
  
  
}


#'
#' @param x
#' @param p
#' @param alpha
#' @param include
#' @param plot
#' 
#' @return 
#' 
#' @example 
#' 
#' 
threshold <-function(x, p, alpha = 0.05, include = 0.75, plot = TRUE) {
  neg <- x <= 0
  f <- p < alpha
  rk <- rank(x)
  nneg <- percTrue(x = -x[neg], f = f[neg], perc = include) 
  npos <- percTrue(x[!neg], f[!neg], perc = include) 
  bs <- bisoft(x, nneg = nneg, npos = npos)
  if (plot) {
    sel <- bs !=0 
    col <- rep("gray", length(p))
    col[p < alpha] <- "red"
    col[sel] <- "green"
    ord <- order(x)
    barplot(x[ord], col = col[ord], border = col[ord])
  }
  bs
}


test.threshold <- function() {
  x <- rnorm(100)
  p <- runif(100)
  
  tx <- threshold(x, p, alpha = 1, include = 0)
  all.equal(x, tx)
  x - tx
  
}


#' 
#' @example 
#' 
sconcord <- function(X, Y, nf=2, dmod = 1, alpha = 0.05, include = 0.8) {
  
  # number of matrices in X, numer of rows in each, numer of columns
  nmat <- length(X) 
  if (is.null(names(X))) 
    names(X) <- paste0("X", 1:nmat)
  
  nr <- sapply(X, nrow)
  nc <- unique(sapply(X, ncol))
  if (length(nc) > 1)
    stop("Number of columns in X need to be the same.")
  
  # sample and feature index
  i.sample <- rep(names(X), each = nc)
  i.feature <- rep(names(X), nr)
  
  # normalization of matrices
  Ynorm <- scale(t(Y), center = TRUE, scale = TRUE)
  Xnorm <- lapply(X, function(x) scale(t(x), center = TRUE, scale = TRUE))
  Xcat <- do.call("cbind", Xnorm)
  
  # define variable
  ys <- gls <- bls <- yloading <- xloading <- xvar <- var <- c()
  
  for (f in 1:nf) {
    cat(paste0("calculating component ", f, "."))
    S <- t(Ynorm) %*% Xcat
    decom <- propack.svd(S, neig = 1, opts = list(kmax = 20))
    load <- decom$v[, 1]
    # dcompLoading <- decom$u[, 1]
    
    ## permutation start
    cat(" Doing permutations .. \n")
    perm <- replicate(30, simplify = FALSE, {
      Sperm <- t(Ynorm[sample(1:nrow(Ynorm)), ]) %*% Xcat
      decomperm <- propack.svd(Sperm, neig = 1, opts = list(kmax = 20))
      list(dcomp = decomperm, S = Sperm)
    })

    loadperm <- sapply(perm, function(x) x$dcomp$v[, 1])
    pm <- pmin(pnorm(q = abs(load), mean = rowMeans(loadperm), sd = rowMads(loadperm), lower.tail = FALSE)*2, 1)
    xloadingSparse <- rep(NA, length(load))
    # loadpermsoft <- loadperm
    for (i in names(X)) {
      ii <- i.sample == i
      xloadingSparse[ii] <- threshold(x = load[ii], p = pm[ii], alpha = alpha, include = include, plot = FALSE)

      # npos <- sum(xloadingSparse[ii] > 0)
      # nneg <- sum(xloadingSparse[ii] < 0)
      # loadpermsoft[ii, ] <- apply(loadperm[ii, ], 2, bisoft, npos = npos, nneg = nneg)
    }

    yloadingcalc <- S %*% xloadingSparse
    # yloadingPerm <- sapply(1:length(perm), function(i) perm[[i]]$S %*% loadpermsoft[, i])
    # py <- pmin(pnorm(q = abs(yloadingcalc),
    #                  mean = rowMeans(yloadingPerm),
    #                  sd = rowMads(yloadingPerm),
    #                  lower.tail = FALSE)*2, 1)
    # yloadingSparse <- threshold(x = yloadingcalc, p = py, alpha = alpha, include = include, plot = TRUE)

    # xloadingSparse <- norm(xloadingSparse)
    # yloadingSparse <- norm(yloadingSparse)
    
    yloadingSparse <- yloadingcalc
    dcompLoading <- norm(yloadingSparse)
    
    
    # xloadingSparse <- decom$v[, 1]
    # yloadingSparse <- decom$u[, 1]
    
    ## permutation end
    xa <- Xcat %*% xloadingSparse
    yb <- Ynorm %*% yloadingSparse
    
    var <- c(var, decom$d[1]^2) ## need to be changed
    ys <- cbind(ys, yb)
    gls <- cbind(gls, xa)
    yloading <- cbind(yloading, yloadingSparse)
    xloading <- cbind(xloading, xloadingSparse)
    
    xai <- lapply(names(X), function(x) {
      ii <- i.feature == x
      Xcat[, ii] %*% matrix(xloadingSparse[ii], ncol = 1)
    })
    
    xai.var <- sapply(xai, crossprod)
    xvar <- cbind(xvar, xai.var)
    bls <- cbind(bls, unlist(xai))
    
    # deflation of S, crossprod matrix, rewrite do not use for loop
    # the classical concordance approach
    # S <- S - tcrossprod(decom$u[, 1]) %*% S
    # # or the same 
    Ynorm <- Ynorm - Ynorm %*% tcrossprod(dcompLoading)
    
  }
  # 
  list(i.sample = i.sample,
       i.feature = i.feature,
       ys = ys,
       ys.norm = apply(ys, 2, norm),
       yloading = yloading,
       gls = gls,
       gls.norm = apply(gls, 2, norm),
       bls = bls,
       bls.norm = apply(bls, 2, tnorm, f=i.sample),
       xloading = xloading,
       xvar = xvar,
       var = var)
}

