library(gdata)
library(shiny)
library(parallel)
library(matrixStats)
library(glmnet)
library(data.table)
library(omic3plus)
library(impute)
library(WebGestaltR)

dat <- readRDS("Res/20170904_concordanceData/nci60_dtp.RDS")
X <- lapply(dat$x, function(x) x$imputed)

vp <- concord(X, dat$y$data, kx = "all", ky = "all", option = "nk", ncomp = 20,
              dmod = 1, center = TRUE, scale = TRUE, pos = FALSE)

loadingy <- data.frame(vp$loading.y)
loadingy <- cbind(drug = rownames(dat$y$data), loadingy)

annotateConcord <- function(obj) {
  
  ag <- sapply(strsplit(rownames(obj$loading.x), " |;"), "[", 2)
  aguniq <- unique(na.omit(ag))
  
  lapply(1:obj$call$ncomp, function(i) {
    cat("================================================================\n\n\n")
    cat(i)
    cat("\n\n\n")
    cat("================================================================\n")
    gs <- data.frame(gene = ag, score = obj$loading.x[, i])
    gs <- gs[!is.na(gs$gene), ]
    gs <- gs[order(abs(gs$score), decreasing = TRUE), ]
    gs <- gs[!duplicated(gs$gene), ]
    
    WebGestaltR(enrichMethod = "GSEA", 
                organism = "hsapiens",
                enrichDatabase="geneontology_Biological_Process_noRedundant",
                minNum=5, maxNum=500,
                interestGene = gs,
                interestGeneType = "genesymbol",
                referenceGene = aguniq, 
                referenceGeneType = "genesymbol",
                collapseMethod = "max",
                is.output = FALSE,
                fdrThr = 0.5)
  })
  
}
annot <- annotateConcord(vp)


i <- 1
plot(vp$score.x[, i], vp$score.y[, i])
abline(a = 0, b = 1)
a <- annot[[i]]


res <- list(concord = vp, annot = annot)
saveRDS(res, file = "Res/20170906_runConcordNCI60/nci60_concord_nosparse_noweight.RDS")
