library(fastmatch)
library(pheatmap)
library(made4)
library(reshape2)
library(fgsea)

prot.tryp <- readRDS(file = "Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
site.tryp <- readRDS(file = "Res/20170110_mq15processLOD2/sites.trypsin.RDS")

prot.tryp.mat <- prot.tryp$fun_prepElnet(prot.tryp, celllines = "all")
site.tryp.mat <- site.tryp$fun_prepElnet(site.tryp, celllines = "all", panel = "NCI60")

an1 <- prot.tryp$annot[rownames(prot.tryp.mat$imputed), ]
an2 <- site.tryp$annot[rownames(site.tryp.mat$imputed), ]
i <- fmatch(an2$`Gene name`, an1$`Gene name`)
ss <- sapply(1:nrow(an2), function(x) {
  cor(prot.tryp.mat$imputed[i[x], ], site.tryp.mat$imputed[x, ])
})
hist(ss, breaks = 100)
names(ss) <- rownames(an2)
ss <- na.omit(ss)
an3 <- an2[names(ss), ]

# function
toGSList <- function(l) {
  ll <- strsplit(l, ";")
  al <- unique(unlist(ll))
  ml <- melt(ll)
  matl <- acast(ml, L1 ~ value, fun.aggregate = length)
  tl <- lapply(1:ncol(matl), function(x) {
    rownames(matl)[matl[, x] == 1]
  })
  names(tl) <- colnames(matl)
  tls
}

dim(an3)
length(ss)

l1 <- c(toGSList(structure(an3$`KEGG name`, names = rownames(an3))),
        toGSList(structure(an3$`GOBP slim name`, names = rownames(an3))),
        toGSList(structure(an3$`GOMF name`, names = rownames(an3))),
        toGSList(structure(an3$`PhosphoSitePlus kinase`, names = rownames(an3))))

l1i <- lapply(l1, function(x) {
  sapply(strsplit(x, " "), "[", 2)
})
l1 <- l1[sapply(l1i, length) >= 5]


l1res <- fgsea(pathways = l1, stats = ss, nperm = 100000, minSize = 5, maxSize = 500, nproc = 20)

l1res$database <- rep(NA, nrow(l1res))
l1res$col <- rep(NA, nrow(l1res))

terms <- unique(unlist(strsplit(an3$`KEGG name`, ";")))
ip <- l1res$pathway %in% terms
l1res$database[ip] <- "KEGG"
l1res$col[ip] <- "orange"

terms <- unique(unlist(strsplit(an3$`GOBP slim name`, ";")))
ip <- l1res$pathway %in% terms
l1res$database[ip] <- "GOBP"
l1res$col[ip] <- "blue"

terms <- unique(unlist(strsplit(an3$`GOMF name`, ";")))
ip <- l1res$pathway %in% terms
l1res$database[ip] <- "GOMF"
l1res$col[ip] <- "cyan"

terms <- unique(unlist(strsplit(an3$`PhosphoSitePlus kinase`, ";")))
ip <- l1res$pathway %in% terms
l1res$database[ip] <- "kinase"
l1res$col[ip] <- "red"

table(l1res$database)
table(is.na(l1res$database))
##


l2i <- lapply(l1res$leadingEdge, function(x) {
  sapply(strsplit(x, " "), "[", 2)
})
l1res


l1res <- l1res[order(l1res$ES), ]



# 
transCol <- function(x) {
  num <- col2rgb(x)
  rgb(num[1, ], num[2, ], num[3, ], alpha = 100, maxColorValue = 255)
}

#
png("Res/__00_figures/site.prot.cor.ES.png", width = 18, height = 12, units = "cm", res = 300)
plot(l1res$ES, pch = 19, col = transCol(l1res$col), cex = -log10(l1res$padj)+0.5, ylab = "Enrichment score")
legend("bottomright", legend = c("KEGG", "GO BP", "GO MF", "Kinase"), col = c("orange", "blue", "cyan", "red"), 
       pch = 19, pt.cex = 2, bty = "n")
dev.off()

png("Res/__00_figures/site.prot.cor.png", width = 18, height = 12, units = "cm", res = 300)
plot(sort(ss), ylab = "Correlation coefficient")
dev.off()

#
plotEnrichment(pathway = l1[["store-operated calcium channel activity"]], stats = ss)
plotEnrichment(pathway = l1[["UL97"]], stats = ss)
plotEnrichment(pathway = l1[["C-X-C chemokine binding"]], stats = ss)
plotEnrichment(pathway = l1[["chemokine binding"]], stats = ss)


pt <- l1res[order(l1res$ES, decreasing = TRUE), ]
