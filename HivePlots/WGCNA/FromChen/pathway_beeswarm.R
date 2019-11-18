outpath <- "/media/kusterlab/users_files/Martin Frejno/postdoc/projects/phosphoproject/res/20180315_wgcna_allsites_allproteins_incnas_nci60/"
eig <- readRDS(file.path(outpath, 'RDS', 'mergedMEs.rds'))$NoImputation

library(beeswarm)


x <- eig$MEantiquewhite4
f <- tryp$nci60_xref$too.short == "LE"

ff <- function(x, f) {
  ll <- split(x, f)
  bxplot(ll, col = "gray", ylab = "module score", frame.plot = FALSE)
  beeswarm(ll, add = TRUE, pch = 19, col = "black", corral = "wrap")
}

tiff("Res/20180316_wgcnaFigures/beeswarm_LE_cAMP.tiff", width = 8, height = 10, units = "cm", res = 300)
ff(eig$MEantiquewhite4, tryp$nci60_xref$too.short == "LE")
dev.off()

tiff("Res/20180316_wgcnaFigures/beeswarm_LE_cellcycle.tiff", width = 8, height = 10, units = "cm", res = 300)
ff(eig$MEsalmon, tryp$nci60_xref$too.short == "LE")
dev.off()

tiff("Res/20180316_wgcnaFigures/beeswarm_RE_hippo.tiff", width = 8, height = 10, units = "cm", res = 300)
ff(eig$MEfirebrick3, tryp$nci60_xref$too.short == "RE")
dev.off()





