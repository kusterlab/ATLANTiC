library(survival)
library(survC1)
library(survHD)
library(survcomp)


fkm <- function(v, s,  ...) {
  srv <- s
  q100 <- seq(0.1, 0.9, by = 0.05)
  ms <- sapply(q100, function(q) {
    md <- quantile(v, na.rm = TRUE, q)
    strata <- rep("low", length(v))
    strata[v > md] <- "high"
    strata <- as.factor(strata)
    if (length(levels(strata)) > 1) {
      mod <- coxph(srv ~ strata)
      pv <- summary(mod)$waldtest["pvalue"]
    } else 
      pv <- 1
    pv
  })
  cp <- q100[which.min(ms)]
  strata <- factor(v > quantile(v, na.rm = TRUE, cp), 
                   levels = c("FALSE", "TRUE"), labels = c("low", "high"))
  plotKM(srv, strata, mark.time=TRUE, ...)
}


cli <- readRDS(file = "Res/20170707_pGR/pgrstainReady.RDS")

f0 <- function(srv) {
  v <- coxph(srv)
  v <- summary(v)
  v$coefficients[, "Pr(>|z|)"][[1]]
}

surv.os <- Surv(cli$`OS-months`, cli$deceased)
coxph(surv.os ~ cli$sperc.mean.phos)
coxph(surv.os ~ cli$sperc.mean.protein)

a <-  cli$sperc.mean.protein + cli$sperc.mean.phos 
coxph(surv.os ~ a)
fkm(v = a, s = surv.os, main = "protein (mean H score)")


