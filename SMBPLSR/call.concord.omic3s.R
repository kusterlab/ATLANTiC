sdata <- readRDS("/media/general/projects/NCI_60_phospho/Chen/Shiny/shinyImage0.3.RDS")
attach(sdata)
library(svd)

fp.nci60 <- readRDS("Res/20170110_mq15processLOD2/nci60.proteinGroups.RDS")
try.nci60 <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
glu.nci60 <- readRDS("Res/20170110_mq15processLOD2/sites.gluC.RDS")

source("R/concordAnalysis/concord.R")
source("R/concordAnalysis/concord.sparse.raw.R")
##

dtp.gi50 <- sdata$drug.nci60$dtp.gi50
any(is.na(dtp.gi50))
dtp.gi50 <- apply(dtp.gi50, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})
dtp.gi50 <- t(dtp.gi50)
any(is.na(dtp.gi50))

prot <- fp.nci60$fun_prepElnet(fp.nci60, celllines = "all")
prot <- prot$imputed
tryp <- try.nci60$fun_prepElnet(try.nci60, celllines = "all", panel = "NCI60", regsitesOnly = FALSE)
tryp <- tryp$imputed
gluc <- glu.nci60$fun_prepElnet(glu.nci60, celllines = 'all', regsitesOnly = FALSE)
gluc <- gluc$imputed


identical(colnames(prot), colnames(tryp))
identical(colnames(prot), paste(colnames(gluc), "NCI60", sep = "_"))
identical(colnames(prot), colnames(dtp.gi50))


nci60 <- list(prot = prot, tryp = tryp, gluc = gluc)


# no sparse
s <- concord(nci60, Y = dtp.gi50, nf = 20)

diag(cor(s$ys, s$gls))

##
xloading <- s$loading
rownames(xloading) <- unlist(lapply(nci60, rownames))
yloading <- s$yloading
rownames(yloading) <- rownames(dtp.gi50)


#
barplot(s$var, ylab = "variance", xlab = "component")

ax <- 1:2
plot(s$ys[, ax[1]], s$gls[, ax[1]], col= sdata$tryp$nci60_xref$color)
text(s$ys[, ax[1]], s$gls[, ax[1]], labels = rownames(s$ys))

plot(s$ys[, ax[2]], s$gls[, ax[2]], col= sdata$tryp$nci60_xref$color)

rownames(s$ys)



# the cell line space
ax <- 1:2
plot(rbind(s$ys.norm[, ax], s$gls.norm[, ax]), col = "white")
points(s$ys.norm[, ax], col = "gray", pch = 20)
arrows(s$ys.norm[, ax[1]], s$ys.norm[, ax[2]], 
       s$gls.norm[, ax[1]], s$gls.norm[, ax[2]], col = "gray", pch = 20)


# text(s$ys[, ax], labels = colnames(gluc), cex = 0.8)

# 
library(ggplot2)
library(gplots)


tiff("Res/20170716_concordFigures/var.tiff", width = 18, height = 12, unit = "cm", res = 300)
barplot(s$var, ylab = "variance", xlab = "component")
dev.off()

qplot(s$var, ylab = "variance", xlab = "component", geom="bar")

dd <- data.frame(cp1 = xloading[, ax[1]], 
                 cp2 = xloading[, ax[2]], 
                 feature = s$i.feature)
tiff("Res/20170716_concordFigures/yloading.tiff", width = 12, height = 12, unit = "cm", res = 300)
qplot(yloading[, ax[1]], yloading[, ax[2]], xlab = "component 1", ylab = "component 2")
dev.off()

tiff("Res/20170716_concordFigures/xloading.tiff", width = 36, height = 12, unit = "cm", res = 300)
qplot(cp1, cp2, data = dd, xlab = "component 1", ylab = "component 2", facets = ~ feature, alpha = 0.5, color = "gray")
dev.off()


ssl <- sort(yloading[, 1], decreasing = TRUE)
barplot(ssl)
ssl[1:6]

rev(ssl)[1:5]


# 
s001 <- xloading[, 1] < -0.01
s001 <- rownames(xloading)[s001]
s001 <- split(s001, substr(names(s001), 1, 1))
s001 <- lapply(s001, function(x) sapply(strsplit(x, " "), "[", 2))
s001 <- lapply(s001, na.omit)

library(VennDiagram)
gh <- venn.diagram(s001, filename = "Res/20170716_concordFigures/venn.tiff", 
                   height = 6, width = 6, units = "cm", res = 300)


intersect(intersect(s001$t, s001$p), s001$g)
intersect(s001$g, s001$p)
setdiff(intersect(s001$t, s001$g), intersect(intersect(s001$t, s001$p), s001$g))


##########################


ssl <- sort(yloading[, 2], decreasing = FALSE)
barplot(ssl)
ssl[1:10]


# 
s001 <- rownames(xloading)[xloading[, 2] < -0.01]
s001 <- split(s001, substr(names(s001), 1, 1))
s001 <- lapply(s001, function(x) sapply(strsplit(x, " "), "[", 2))
s001 <- lapply(s001, na.omit)


intersect(intersect(s001$t, s001$p), s001$g)
intersect(s001$g, s001$p)
setdiff(intersect(s001$t, s001$g), intersect(intersect(s001$t, s001$p), s001$g))

library(VennDiagram)
gh <- venn.diagram(s001, filename = "Res/20170716_concordFigures/venn_p2.tiff", 
                   height = 6, width = 6, units = "cm", res = 300)


write.table(unique(t(t(unlist(s001)))), file = "temp.txt", col.names = FALSE, 
            row.names = FALSE, quote = FALSE, sep ="\t")

# plot(xloading[, ax])
ax1 <- 2
par(mar = c(10, 4, 0.1, 0.1))
barplot(sort(xloading[, ax1], decreasing = FALSE)[1:30], las = 2)
barplot(sort(yloading[, ax1], decreasing = FALSE)[1:10], las = 2)


for (i in 1:20) {
  print(i)
  xx <- yloading[, i]
  print(sort(xx)[1:3])
  print(sort(xx, decreasing = TRUE)[1:3])
  print("============")
  
}

ss <- yloading[, 10]
tiff("bar.yloading.10.tiff", width = 18, height = 12, unit = "cm", res = 300)
barplot(sort(ss), ylab = "component 10")
dev.off()


barplot(sort(xloading[, 10]))

s001 <- rownames(xloading)[xloading[, 10] < -0.01]
s001 <- split(s001, substr(names(s001), 1, 1))
s001 <- lapply(s001, function(x) sapply(strsplit(x, " "), "[", 2))
s001 <- lapply(s001, na.omit)

gh <- venn.diagram(s001, filename = "Res/20170716_concordFigures/venn_p10.tiff", 
                   height = 6, width = 6, units = "cm", res = 300)


intersect(intersect(s001$t, s001$p), s001$g)
intersect(s001$g, s001$p)
setdiff(intersect(s001$t, s001$g), intersect(intersect(s001$t, s001$p), s001$g))

# names(sort(xloading[, ax1], decreasing = FALSE)[1:400])
# 
# ##
# library(rbokeh)
# library(scatterD3)
# scatterD3(yloading[, ax[1]], yloading[, ax[2]], point_opacity = 0.5, tooltips = rownames(yloading), hover_size = 2)

##
sr <- sconcord(nci60, Y = dtp.gi50, nf = 5, alpha = 0.05, include = 0.75)




