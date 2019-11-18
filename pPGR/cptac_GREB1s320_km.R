library(survival)
library(survC1)
library(survHD)
library(survcomp)

clinical <- readRDS("Res/20180411_breastCptacSurv/clin.RDS")
drug <- readRDS("Res/20180411_breastCptacSurv/clin.drug.RDS")
  
hist(clinical$GREB1S320)

srv <- Surv(time = as.numeric(clinical$OS), event = c(clinical$vital_status))
plotKM(y = srv, strata = as.factor(clinical$GREB1S320 > 2), mark.time=TRUE, main = "OS")
