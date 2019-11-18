require(data.table)
require(reshape2)
require(drc)
source("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/pubFigures2/R/colorScheme.R")
source('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/src/plotFACI_hack.R')
# source('/media/msdata5/users_files/martin/phd/data/ms/Drug_sensitivity/synergy/chou8491.R')

folder <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/followups/met_ron/'
files <- list.files(folder, full.names = T, pattern = '.txt')
datlist <- lapply(files, function(i) fread(i, stringsAsFactors = F, integer64 = 'double'))
names(datlist) <- unlist(lapply(datlist, function(i) colnames(i)[1]))
concs <- names(datlist)
celllines <- unique(unlist(lapply(datlist, function(i) colnames(i))))[!unique(unlist(lapply(datlist, function(i) colnames(i))))%in%concs]
lapply(concs, function(i) datlist[[i]][,c(i):=10^eval(as.name(i))])
lapply(concs, function(i) setkeyv(datlist[[i]], i))
# meanconc <- sapply(concs, function(i) datlist[[i]][,mean(eval(as.name(i)),na.rm = T)])
# concs <- names(meanconc)
# concs <- c(concs[grep('mix', ignore.case = T, concs)],names(sort(meanconc[grep('mix', ignore.case = T, names(meanconc),invert = T)], decreasing = T)))
# datlist <- datlist[concs]
singles <- grep('mix', x = concs, ignore.case = T, files, invert = T, value = T)
singidx <- grep('mix', x = concs, ignore.case = T, files, invert = T)
conctab <- data.table(sapply(concs, function(i) datlist[[i]][,eval(as.name(i))]))
conctab[,colnames(conctab)[!colnames(conctab)%in%singles]:=eval(as.name(singles[1]))+eval(as.name(singles[2]))]
datlist[[concs[!concs%in%singles]]][,concs[!concs%in%singles]:=conctab[,eval(as.name(concs[!concs%in%singles]))]]
for (i in seq_along(concs)) {
  colnames(datlist[[i]])[!colnames(datlist[[i]])%in%c(concs, 'comb')] <- paste(seq(1,length(colnames(datlist[[i]])[!colnames(datlist[[i]])%in%c(concs, 'comb')]),1), colnames(datlist[[i]])[!colnames(datlist[[i]])%in%c(concs, 'comb')], sep = '_')
}
datlist_l <- lapply(seq_along(concs), function(i) melt(datlist[[i]], id.vars = concs[i])[,variable:=gsub('[0-9]_', '', variable)])
names(datlist_l) <- concs
lapply(concs, function(i) setkeyv(datlist_l[[i]], c('variable', i)))
concdt <- data.table(sapply(concs, function(i) datlist_l[[i]][,eval(as.name(i))]))
props <- sapply(singles, function(i) concdt[,mean(eval(as.name(i))/eval(as.name(concs[!concs%in%singles])), na.rm=T)])
# ratio <- ifelse(signif(mean(concdt[,V1/V2]), 5)<1,signif(mean(concdt[,V1/V2]), 5),1/signif(mean(concdt[,V1/V2]), 5))
modellist <- lapply(celllines, function(j) lapply(concs, function(i) try(drm(formula = value ~ eval(as.name(i)), data = datlist_l[[i]][variable==j,], curveid = as.factor(variable), fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50"))), silent = T)))
names(modellist) <- celllines
for (i in seq_along(modellist)) {
  names(modellist[[i]]) <- concs
}
png(filename = file.path(folder, 'CIplot.png'), width = 23.7, height = 12.5, units = "cm", pointsize = 10, res = 300)
par(mfrow=c(1,length(celllines)))
lapply(celllines, function(j) plotFACI(effList = CIcompX(mixProp = props[1], modelList = modellist[[j]][concs], EDvec = seq(1,99,1), EDonly = F), indAxis = 'ED', caRef = T, showPoints = F, ylim = c(0,3), main = j))
dev.off()
lapply(celllines, function(j) CIcomp(mixProp = props[1], modelList = modellist[[j]][concs], EDvec = seq(0,100,1)))

setEPS()
postscript(file.path(folder, 'CIplot.eps'), width = 5.8 /2.539998, height = 5.8/2.539998, pointsize = 10, family = 'TUM Neue Helvetica 55 Regular', encoding = 'Greek')
# par(mfrow=c(1,length(celllines)))
par(mar = c(2.5, 3.1, 2.3, 1.2))
plotFACI_FIX(effList = CIcompX(mixProp = props[1], modelList = modellist[['HDC8']][concs], EDvec = seq(1,100,2), EDonly = F), indAxis = 'ED', caRef = T, showPoints = F, ylim = c(0,3.1), main = '', font.main = 1, cex.main = 1, xaxt = 'n', yaxt = 'n', col = 'white', bty = 'n', xlab = NA, ylab = NA)
rect(xleft = -10, ybottom = -1, xright = 110, ytop = 1, density = 100, col = al2hex(.cs$tColor(x = .cs$colorScheme(scheme = 'c3')[1], alpha = 50)), border = al2hex(.cs$tColor(x = .cs$colorScheme(scheme = 'c3')[1], alpha = 50)))
text(x = 50, y = 0.3, adj = c(0.5,0.5), labels = 'Synergy', col = .cs$colorScheme(scheme = 'c3')[1])
rect(xleft = -10, ybottom = 1, xright = 110, ytop = 4, density = 100, col = al2hex(.cs$tColor(x = .cs$colorScheme(scheme = 'c3')[2], alpha = 50)), border = al2hex(.cs$tColor(x = .cs$colorScheme(scheme = 'c3')[2], alpha = 50)))
text(x = 50, y = 2, adj = c(0.5,0.5), labels = 'Antagonism', col = .cs$colorScheme(scheme = 'c3')[2])
plotFACI_FIX(effList = CIcompX(mixProp = props[1], modelList = modellist[['HDC8']][concs], EDvec = seq(1,100,2), EDonly = F), indAxis = 'ED', caRef = T, showPoints = F, ylim = c(0,3.1), main = '', font.main = 1, cex.main = 1, xaxt = 'n', yaxt = 'n', add = T, bty = 'n')
box(lwd = 1)
axis(1, at = seq(0,100,by= 40), line = 0, cex.axis = 1, padj = 0.5, tck=-0.02, mgp = c(3,0,0))
axis(1, at = seq(20,100,by= 40), line = 0, cex.axis = 1, padj = 0.5, tck=-0.02, mgp = c(3,0,0))
axis(2, line = 0, cex.axis = 1, padj = 0.5, tck=-0.02, mgp = c(3,0.5,0), las = 2)
abline(h = 1, lty = 2, lwd = 1, col = .cs$colorScheme(scheme = 'c3')[3])
legend('topright', legend = 'Additivity', lty = 2, col = .cs$colorScheme(scheme = 'c3')[3], bty = 'n', text.col = .cs$colorScheme(scheme = 'c3')[3])
mtext(side = 1, text = "Fraction affected", cex = 1, line = 1.1)
mtext(side = 2, text = "Combination index", cex = 1, line = 2)
mtext(side = 3, text = "Combination treatment\ntargeting MET and MST1R - CI", cex = 1, line = 0.2, adj = 0.5)
dev.off()

CIs <- CIcompX(mixProp = props[1], modelList = modellist[['HDC8']][concs], EDvec = seq(0,100,2), EDonly = F)[['CAx']]
EDs <- CIcompX(mixProp = props[1], modelList = modellist[['HDC8']][concs], EDvec = seq(0,100,2), EDonly = F)[['EDvec']]

lines(100-EDs, CIs[,"combInd"], col = 2)

CIs[which.min(CIs[,"CAdiffp"]),]
EDs[which.min(CIs[,"CAdiffp"])]
CIs[which(EDs==50),]
# 
# 
# plot(dat$cor, y = -log10(dat$pval), col = "white", cex = 1, axes = FALSE, pch = 19, xlab = "", ylab = "") # , ylim = c(0,max(-log10(dat$pval)))
# # abline(h = 1:11, col = "gray", lty = 3)
# box()
# points(dat$cor[dat$cor>0], y = -log10(dat$pval)[dat$cor>0], col = col['TUMred'], cex = 1.1, pch = 19)
# points(dat$cor[dat$cor<0], y = -log10(dat$pval)[dat$cor<0], col = col['TUMblue'], cex = 1.1, pch = 19)
# points(dat$cor, y = -log10(dat$pval), col = "white", cex = 0.9, pch = "_")
# # mtext(side = 2, text = 1:11, line = 0.2, at = 1:11, las = 2, cex = 0.6)
# 


# combInd <- (mixProp * eseVec[1]/eseVec[2]) + ((1 - mixProp) * eseVec[1]/eseVec[3])
# extCI <- function()

# lapply(datlist, function(i) i[,singles[!singles %in% colnames(i)]:=do.call(cbind, lapply(singles[!singles %in% colnames(i)], function(j) datlist[[which(concs%in%j)]][,list(eval(as.name(j)))]))])
# lapply(singles, function(i) datlist[[which(concs%in%i)]][,singles[!singles%in%i]:=0])
# # lapply(concs, function(i) datlist[[which(concs%in%i)]][,comb:=i])
# datlist[[which(!concs%in%singles)]][,concs[!concs%in%singles]:=NULL]
# dat <- data.frame(rbindlist(datlist, use.names = T), check.names = F)
# colnames(dat)[!colnames(dat)%in%c(singles, 'comb')] <- paste(seq(1,length(colnames(dat)[!colnames(dat)%in%c(singles, 'comb')]),1), colnames(dat)[!colnames(dat)%in%c(singles, 'comb')], sep = '_')
# # dat_dt <- data.table(dat)
# # dat_dt[,paste0('10^', singles[1]):=10^eval(as.name(singles[1]))]
# # dat_dt[,paste0('10^', singles[2]):=10^eval(as.name(singles[2]))]
# dt <- melt(data.table(dat), id.vars = c(singles,'comb'))
# dt[,variable:=gsub('[0-9]*_', '', variable)]
# dt[,variable:=paste(comb,variable,sep = '_')]
# dt[,comb:=NULL]
# 
# # m <- dat[eval(as.name(singles[1]))!=0 | eval(as.name(singles[2]))!=0,unique(10^eval(as.name(singles[2]))/10^eval(as.name(singles[1])))]
# # m <- unique(syn1$UNC569/syn1$RDEA119)[4]
# # drug.combination.chou(data.frame(dat[,c(singles,"value"),with=F]), d2.d1 = 0.01)
