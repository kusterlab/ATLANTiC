library(beeswarm)
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
}


cli <- readRDS(file = "Res/20170707_pGR/pgrstainReady.RDS")
plot(cli$h.score.mean.phos, cli$h.score.mean.protein)
cor(cli$h.score.mean.phos, cli$h.score.mean.protein, use = "pair")
table(cli$PR.STATUS)
beeswarm(cli$h.score.mean.protein ~ cli$PR.STATUS, corral  = "wrap")
abline(h = 85)


cli <- cli[which(cli$h.score.mean.protein > 85 & cli$PR.STATUS == 1), ]
ss <- cli$`OS-months`
ds <- cli$deceased
ss <- Surv(time = ss, event = ds)

fkm(v = cli$h.score.mean.protein, s = ss, main = "protein (mean H score)")
fkm(v = cli$h.score.mean.phos, s = ss, main = "phospho-site (mean H score)")

fkm(cli$h.score.max.protein, s = ss, main = "protein (max H score)")
fkm(cli$h.score.max.phos, s = ss, main = "phospho-site (max H score)")

fkm(cli$sintens.mean.protein, s = ss, main = "protein (mean intensity)")
fkm(cli$sintens.mean.phos, s = ss, main = "phospho-site (mean intensity)")

fkm(cli$sintens.max.protein, s = ss, main = "protein (max intensity)")
fkm(cli$sintens.max.phos, s = ss, main = "phospho-site (max intensity)")

fkm(cli$sperc.mean.protein, s = ss, main = "protein (mean perctage)")
fkm(cli$sperc.mean.phos, s = ss, main = "phospho-site (mean percetage)")

fkm(cli$sperc.max.protein, s = ss, main = "protein (max percetage)")
fkm(cli$sperc.max.phos, s = ss, main = "phospho-site (max percetage)")

# dev.off()










