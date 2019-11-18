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


# survtime <- dat$OS
# evt <- dat$`OS-Status`

srv <- Surv(time = as.numeric(survtime), event = evt)
plotKM(y = srv, strata = as.factor(dat$`SAMHD1.Score.(binary)`), mark.time=TRUE, main = "SAMHD1 (OS)")

srv <- Surv(time = as.numeric(survtime), event = evt)
plotKM(y = srv, strata = as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`), 
       mark.time=TRUE, main = "AK1 (OS)")


i <- dat$`EFS-Status.censored.for.SCT` %in%  c("event (relapse)", "event (relapse) ", "no event", "no-event")
table(i)
srv <- Surv(time = as.numeric(survtime)[i], event = evt[i])
plotKM(y = srv, strata = as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`[i]), 
       mark.time=TRUE, main = "AK1 (SAMHD1 high only [OS])")
# dev.off()



