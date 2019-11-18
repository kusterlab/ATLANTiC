library(survival)
library(survC1)
library(survHD)
library(survcomp)

cli <- readRDS(file = "Res/20170707_pGR/pgrstainReady.RDS")


v <- coxph(Surv(cli$`OS-months`, event = cli$deceased) ~ cli[[20]] + cli$HER2.STATUS)
v$coefficients


f0 <- function(srv) {
  v <- coxph(srv)
  v <- summary(v)
  v$coefficients[, "Pr(>|z|)"][[1]]
  # v$logtest[["pvalue"]]
}

coxtest <- function(x, s) {
  res <- c(
    x = f0(Surv(s, event = cli$deceased) ~ x),
    age = f0(Surv(s, event = cli$deceased) ~ x + cli$age),
    pT = f0(Surv(s, event = cli$deceased) ~ x + cli$pT),
    pN = f0(Surv(s, event = cli$deceased) ~ x + cli$pN),
    pN_status = f0(Surv(s, event = cli$deceased) ~ x + cli$pN_status),
    grade = f0(Surv(s, event = cli$deceased) ~ x + cli$grading),
    
    er_status = f0(Surv(s, event = cli$deceased) ~ x + cli$ER.STATUS),
    pr_status = f0(Surv(s, event = cli$deceased) ~ x + cli$PR.STATUS),
    her2_status = f0(Surv(s, event = cli$deceased) ~ x + cli$HER2.STATUS)
  )
  res
}

dfs <- sapply(cli[20:31], function(x) {
  coxtest(x = x, s = cli$`DFS-months`)
})

os <- sapply(cli[20:31], function(x) {
  coxtest(x = x, s = cli$`OS-months`)
})

library(pheatmap)
pheatmap(-log10(dfs), cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, filename = "Res/20170707_pGR/DFS.coxph.png")
pheatmap(-log10(os), cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, filename = "Res/20170707_pGR/OS.coxph.png")
