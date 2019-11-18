library(openxlsx)
library(survival)
library(survC1)
library(survHD)
library(survcomp)

dat <- read.xlsx("Dat/AK1/AML_AK1 IHC_TO.xlsx")

table(dat$`AK1.Score.clone.OAAB17548.1:50`)
table(dat$`AK.Score.(binary).clone.OAAB17548.1:50`)
############################################################

v <- coxph(srv ~  as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`)+ grepl("46", dat$`-ryotype`))
v$coefficients
exp(v$coefficients)

v <- summary(v)

f0 <- function(srv) {
  v <- coxph(srv)
  v <- summary(v)
  x <- v$coefficients[, "Pr(>|z|)"]
  x <- c(x[1], min(x[-1]))
  names(x) <- c("AK1low", "fac")
  x[is.infinite(x)] <- NA
  x
}

cf <- function(srv) {
  cbind(
    AK1only = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`)),
    age = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + dat$Age),
    sex = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + dat$Sex),
    Blasts.BM.at.diagnosis = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + dat$`Blasts.%.BM.at.diagnosis`),
    Therapy.Induction.2.Consolidation = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + dat$`Therapy:.Induction.2./.Consolidation`),
    Status.CR.after.induction.therapy = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + c(dat$Status.CR.after.induction.therapy == "CR")),
    Stem.cell.transplantation = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + dat$Stem.cell.transplantation),
    Initial.WBC.nl = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + dat$`Initial.WBC/nl`),
    FAB = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + dat$FAB),
    NPM = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + (dat$NPM == "negative")),
    MLL = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + (dat$MLL == "negative")),
    FLT3 = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + (dat$FLT3 == "wt")),
    inv16 = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + dat$`inv(16)`),
    t.8.21 = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + dat$`t(8;21)`),
    ckit = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + (dat$`c-Kit` == "negative")),
    CEBPalpha = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + dat$CEBP.alpha),
    ELN.cytogenetic.risk = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + tolower(dat$ELN.cytogenetic.risk)),
    karyotype = f0(srv ~ as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`) + grepl("46", dat$`-ryotype`))
  )
}


# ==================================== Overall survival ========================================== 

png("Res/20170821_AK1survanalysis/OS.png", width = 28, height = 16, units = "cm", res = 200)
layout(matrix(1:2, 1, 2))
srv <- Surv(time = as.numeric(dat$OS), event = dat$`OS-Status` == "dead")
plotKM(y = srv, strata = as.factor(dat$`AK.Score.(binary).clone.OAAB17548.1:50`), mark.time=TRUE, main = "OS")
# plotKM(y = srv, strata = as.factor(dat$`AK1.Score.clone.OAAB17548.1:50`), mark.time=TRUE)
par(mar = c(15, 4, 4, 1))
val <- -log10(cf(srv))
barplot(val, beside = TRUE, las = 2, ylab = "p-value (coxph model; -log10)", 
        legend.text = c("AK1 binary", "other factors"), args.legend = list(x = "bottomright"))
abline(h = -log10(0.05), lty = 2, col = "red")
abline(h = val[1, 1], lty = 2, col = "green")
dev.off()

