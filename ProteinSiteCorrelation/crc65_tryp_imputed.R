library(parallel)

tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
crc <- readRDS("Res/20170110_mq15processLOD2/crc65.proteinGroups.RDS")

sc <- setdiff(colnames(crc$ibaq), "CoCM-1_CRC65")
expr_site <- tryp$fun_prepElnet(tryp, panel = "CRC65", celllines = sc)$imputed
expr_prot <- crc$fun_prepElnet(crc, celllines = sc)$imputed

annot_site <- tryp$annot[rownames(expr_site), ]
annot_prot <- crc$annot[rownames(expr_prot), ]

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
tiff("../Manuscript/supplementFigures/supp_site_prot_cor_crc65_imputed.tiff", 
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
# negcor <- fres[fres$r < -0.5 & fres$sitemax > 8, ]
# f1 <- function(x) sapply(strsplit(x, ";"), "[", 1)
# gn1 <- f1(unique(annot_prot[fres$prot[fres$r < -0.25], "Gene name"]))
# gn2 <- f1(unique(annot_prot[fres$prot[fres$r > -0.25 & fres$r < 0.25], "Gene name"]))
# gn3 <- f1(unique(annot_prot[fres$prot[fres$r > 0.25], "Gene name"]))
# 
# write.table(gn1, file = "tmp.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
