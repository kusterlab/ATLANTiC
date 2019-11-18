library(mogsa)

X <- list(prot = prot, tryp = tryp, gluc = gluc)
Y <- dmat

col <- fp.nci60$xref$color #

rr <- sconcord(X, Y, nf = 6)

diag(cor(rr$scoreXall, rr$scoreY))
ax <- 5:6
plot(rr$scoreXall[, ax], col = fp.nci60$xref$color, pch = 19, xlim = c(-80, 60))
plot(rr$scoreY[, ax], col = fp.nci60$xref$color, pch = 19, xlim = c(-180, 60))
text(rr$scoreY[, ax], labels = rownames(rr$scoreY))

plot(rr$loadingXs[, 3:4])
plot(rr$loadingY0[, 3:4])

drug.nci60$dtp.annot[rr$loadingY0[, 3] < 0, ]
rownames(prot)[rr$loadingXs[1:nrow(prot), 3] < 0]


# plot(rr[[1]]$t, rr[[1]]$t0)
# rownames(dmat)[rr[[1]]$t0 != 0]
# drug.nci60$dtp.annot[rr[[1]]$t0 != 0, ]
# 
# plot(rr[[1]]$t0[, 1], rr[[2]]$t0[, 1], col = col, pch = 19)
# plot(rr[[3]]$t0[, 1], rr[[4]]$t0[, 1], col = col, pch = 19)
# 
# r <- rr[[5]]$corsum
# fg(rr[[2]]$corsum)
# barplot(rr[[4]]$corsum)
# which.max(rr[[5]]$corsum)
# rr[[1]]$t

