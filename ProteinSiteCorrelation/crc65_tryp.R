library(parallel)
library(psych)

tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
crc <- readRDS("Res/20170110_mq15processLOD2/crc65.proteinGroups.RDS")

sc <- setdiff(colnames(crc$ibaq), "CoCM-1_CRC65")
expr_site <- tryp$crc65_expr[, sc]
expr_prot <- crc$ibaq[, sc]

annot_site <- tryp$annot
annot_prot <- crc$annot

gene_idx_site <- sapply(strsplit(annot_site$`Gene name`, ";"), "[", 1)
gene_idx_prot <- sapply(strsplit(annot_prot$`Gene name`, ";"), "[", 1)

sharedGenes <- intersect(gene_idx_site, gene_idx_prot)

res <- mclapply(sharedGenes, function(x) {
  i_prot <- which(gene_idx_prot == x)
  i_site <- which(gene_idx_site == x )
  if (length(i_prot) < 1 | length(i_site) < 1)
    return(NULL)
  m1 <- expr_site[i_site, , drop = FALSE]
  m2 <- expr_prot[i_prot, , drop = FALSE]
  ret <- corr.test(t(m1), t(m2), use = "pair", ci = FALSE, adjust = "none")
  
  rmax <- rowMaxs(m1, na.rm = TRUE)
  
  df <- data.frame(r = c(ret$r), 
                   n = c(ret$n),
                   t = c(ret$t),
                   p = c(ret$p), 
                   sitemax = rmax, 
                   site = rownames(ret$r), 
                   prot = rep(colnames(ret$r), each = nrow(ret$r)), 
                   stringsAsFactors = FALSE)
  df
}, mc.cores = 20)
res <- do.call(rbind, res)

fres <- res[res$n >= 7, ]
hist(res$r)
hist(res$r[res$n >= 7])


## plot
tiff("../Manuscript/supplementFigures/supp_site_prot_cor_crc65.tiff", 
     width = 8.7, height = 8,7, units = "cm", res = 300)
layout(matrix(1:4, 2, 2), widths = 1:2, heights = 1:2)
plot(x = 0, y = 0, pch = ".", frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "")

par(mar = c(3, 2, 0, 0.5))
px <- hist(fres$sitemax, breaks = 100, plot = FALSE)
barplot(height = -px$counts, horiz = TRUE, axes = FALSE, border = "gray")
axis(side = 1, labels = seq(0, 800, by = 200), at = -seq(0, 800, by = 200))
mtext(side = 1, "Counts", line = 2)

par(mar = c(0.1, 0, 0.5, 3))
px <- hist(fres$r, breaks = 100, plot = FALSE)
barplot(px$counts, axes = FALSE, border = "gray")
axis(4)
mtext(side = 4, "Counts", line = 2)

par(mar = c(3, 0, 0, 3))
# plot(fres$r, fres$sitemax, pch = 20, xlim = c(-1, 1), axes = FALSE, cex = 0.5, 
# col = rgb(190, 190, 190, alpha = 100, maxColorValue = 255))
smoothScatter(fres$r, fres$sitemax, axes = FALSE)
box()
axis(1)
mtext("Max intensity of p-sites", side = 4, line = 2)
axis(4)
mtext("Pearson R", side = 1, line = 2)
dev.off()

## 

# 
# ###
f1 <- function(x) sapply(strsplit(x, ";"), "[", 1)
gn2 <- f1(unique(annot_prot[fres$prot[fres$r < 0.3], "Gene name"]))
gn3 <- f1(unique(annot_prot[fres$prot[fres$r > 0.7], "Gene name"]))

library(WebGestaltR)


em2 <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                   enrichDatabase="geneontology_Biological_Process",
                   interestGene=gn2,
                   interestGeneType="genesymbol", referenceSet = "genome",
                   minNum=10, maxNum=500,
                   fdrMethod="BH", sigMethod="fdr",
                   fdrThr=0.1, topThr=10,
                   dNum=20, perNum=1000, lNum=20,
                   is.output=FALSE)

em3 <- WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                   enrichDatabase="geneontology_Biological_Process",
                   interestGene=gn3,
                   interestGeneType="genesymbol", referenceSet = "genome",
                   minNum=10, maxNum=500,
                   fdrMethod="BH", sigMethod="fdr",
                   fdrThr=0.1, topThr=10,
                   dNum=20, perNum=1000, lNum=20,
                   is.output=FALSE)

write.table(em2[1:10, c(1:2, 4:9)], file = "tmp.txt", 
            col.names = TRUE, row.names = FALSE, quote =FALSE, sep = "\t")
write.table(em3[1:10, c(1:2, 4:9)], file = "tmp2.txt", 
            col.names = TRUE, row.names = FALSE, quote =FALSE, sep = "\t")

