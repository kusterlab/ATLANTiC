library(openxlsx)
library(survival)
library(survC1)
library(survHD)
library(survcomp)

dat <- read.xlsx("Dat/AK1/AML_AK1_SAMHD1_IHC_TO.xlsx")
dat <- dat[dat$`SAMHD1.Score.(binary)` != "no tissue", ] 

table(dat$`AK1.Score.clone.OAAB17548.1:50`)
table(AK1 = dat$`AK.Score.(binary).clone.OAAB17548.1:50`, 
      SAMHD1 = dat$`SAMHD1.Score.(binary)`)

#####
survtime <- dat$OS
evt <- dat$`OS-Status`

dir.create(path = "Res/201805014_AK1SAMHD1/")
png("Res/201805014_AK1SAMHD1/OS.png", width = 26, height = 8, units = "cm", res = 300)
layout(matrix(1:3, 1, 3))
srv <- Surv(time = as.numeric(survtime), event = evt == "dead")
plotKM(y = srv, strata = as.factor(dat$`SAMHD1.Score.(binary)`), mark.time=TRUE, main = "SAMHD1 (OS)")

srv <- Surv(time = as.numeric(survtime), event = evt == "dead")
plotKM(y = srv, strata = as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`), 
       mark.time=TRUE, main = "AK1 (OS)")

i <- dat$`SAMHD1.Score.(binary)` == "high"
table(i)
srv <- Surv(time = as.numeric(survtime)[i], event = evt[i] == "dead")
plotKM(y = srv, strata = as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`[i]), 
       mark.time=TRUE, main = "AK1 (SAMHD1 high only [OS])")
dev.off()

srv <- Surv(time = as.numeric(survtime), event = evt == "dead")
plotKM(y = srv, strata = as.factor(dat$`AK1.Score.clone.OAAB17548.1:50`), mark.time=TRUE, main = "OS")

i <- dat$SAMHD1.Score == "3"
table(i)
srv <- Surv(time = as.numeric(survtime)[i], event = evt[i] == "dead")
plotKM(y = srv, strata = as.factor(dat$`AK1.Score.clone.OAAB17548.1:50`[i]), 
       mark.time=TRUE, main = "OS")


table(samhd1 = dat$SAMHD1.Score, 
      ak1 = dat$`AK1.Score.clone.OAAB17548.1:50`)

srv <- Surv(time = as.numeric(survtime), event = evt == "dead")
plotKM(y = srv, strata = as.factor(dat$`SAMHD1.Score.(binary)`), mark.time=TRUE, main = "OS")

srv <- Surv(time = as.numeric(survtime), event = evt == "dead")
plotKM(y = srv, strata = as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`), mark.time=TRUE, main = "OS")

srv <- Surv(time = as.numeric(survtime), event = evt == "dead")
v <- paste(dat$`AK.Score.(binary).clone.OAAB17548.1:50`, dat$`SAMHD1.Score.(binary)`)
v[ !v %in% names(which(table(v) > 5)) ] <- NA
plotKM(y = srv, strata = as.factor(v), mark.time=TRUE, main = "OS")
