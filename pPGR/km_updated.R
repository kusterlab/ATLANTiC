library(survival)
library(survC1)
library(survHD)
library(survcomp)

fkm <- function(v, s,  ...) {
  srv <- s
  q100 <- seq(0.1, 0.9, by = 0.05)
  ms <- sapply(q100, function(q) {
    md <- quantile(v, na.rm = TRUE, q)
    # strata <- factor(v > md, levels = c("FALSE", "TRUE"), labels = c("low", "high"))
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
  v > quantile(v, na.rm = TRUE, cp)
}

## ----------- use clinical classification ------------------


cli <- readRDS(file = "Res/20170707_pGR/pgrstainReady.RDS")
surv.os <- Surv(cli$`OS-months`, cli$deceased)

tiff("Res/20190109_newPGR/KM_clinicianClass.tiff", height = 20, width = 13, units = "cm", res = 300)
layout(matrix(1:2, 2, 1))
iv <- fkm(v = cli$PR.STATUS, s = surv.os, main = " PGR (Clinician's classification)", censor.at = 120)
cli1 <- cli[cli$PR.STATUS == 1, ]
surv.os1 <- Surv(cli1$`OS-months`, cli1$deceased)
iv <- fkm(v = cli1$h.score.mean.phos, s = surv.os1, main = "PGR pS20 (H score)", censor.at = 120)
dev.off()

