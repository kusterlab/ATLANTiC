library(PTMotif)
library(parallel)
source("pubFigures2/R/colorScheme.R")

nci60 <- readRDS(file = "Res/20180515_kinaseSub/filteredDF.NCI60.RDS")
nci60$panel <- "NCI60"
crc65 <- readRDS(file = "Res/20180515_kinaseSub/filteredDF.CRC65.RDS")
crc65$panel <- "CRC65"
subs <- rbind(nci60, crc65)
write.table(subs, file = "Res/20180515_kinaseSub/supp_kinaseSubstrateTable.txt", 
            col.names = TRUE,row.names = FALSE, quote= FALSE, sep = "\t")
subs <- split(subs, subs$kinase)

###
integ <- readRDS("Res/20160803_kinaseSubstrateDatabaseParse/integDatabase.RDS")
km <- read.delim("Res/20160803_kinaseSubstrateDatabaseParse/kinaseMap.txt", sep = "=", header = FALSE, stringsAsFactors = FALSE)
km$V1 <- trimws(km$V1)
km$V2 <- trimws(km$V2)
intersect(as.character(integ$kinase), km$V1)
integ$kinaseMapTo <- km$V2[match(as.character(integ$kinase), km$V1)]
write.table(integ, file = "Res/20180515_kinaseSub/supp_integKinaseSubDB.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


###
kin <- as.character(integ$kinaseMapTo)
kin[is.na(kin)] <- as.character(integ$kinase)[is.na(kin)]
kin <- data.frame(ki = kin,
                 sub = as.character(integ$sub_geneName), 
                 stringsAsFactors = FALSE
                 )
sl <- strsplit(kin$ki, ";")
kin <- data.frame(ki = unlist(sl), 
                  rep(kin$sub, sapply(sl, length)), 
                  stringsAsFactors = FALSE)
kin <- unique(kin)
integ <- split(kin$rep.kin.sub..sapply.sl..length.., kin$ki)

### actual mapping
ni <- intersect(names(subs), names(integ))
v <- lapply(ni, function(nn) {
  x <- integ[[nn]]
  y <- subs[[nn]]
  i <- as.character(y$substrate) %in% x
  y[i, ]
})
names(v) <- ni
v0 <- lapply(v, function(x) length(unique(as.character(x$substrate))))
vt <- do.call(rbind, v)
write.table(vt, file = "Res/20180515_kinaseSub/supp_overlapTable.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


##
ol <- unlist(v0)
rb <- rbind(overlap = ol,
            corOnly = sapply(subs[ni], nrow) - ol, 
            dbOnly = sapply(integ[ni], length) - ol
            )
write.table(t(rb), file = "Res/20180515_kinaseSub/supp_overlapNumbers.txt", 
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")


rb <- rb[, colSums(rb[c(1, 3), ]) > 50]
ord <- order(rb[1, ], decreasing = TRUE)
rb <- rb[, ord]          

setEPS()
postscript("pubFigure3/figures/newfig2_kinaseSubstrateOverlap.eps", 
           width = 17.4/2.539998, height = 7/2.539998, pointsize = 10,
           family = 'TUM Neue Helvetica 55 Regular', encoding = 'Greek')
par(mar = c(3, 4, 3, 5), lwd = 0.01, xpd = TRUE)
bp <- barplot(sweep(rb, 2, colSums(rb), "/"), beside = FALSE, width = colSums(rb)+50, 
              # border = .cs$colorScheme("c3"), 
              border = "white", 
              las = 2, lwd = 0.1, cex.names = 0.3, 
              col = .cs$colorScheme("c3"))
mtext(rb[1, ], line = 0.5, side = 3, at = bp, cex = 0.3, las = 2)
mtext(rb[2, ], line = 1.0, side = 3, at = bp, cex = 0.3, las = 2)
mtext(rb[3, ], line = 2, side = 3, at = bp, cex = 0.3, las = 2)

legend(x = max(bp)*1.05, y = 1, legend = rownames(rb), col = .cs$colorScheme("c3"), pch = 15, bty = "n")

dev.off()


setEPS()
postscript("pubFigure3/figures/newfig2_kinaseSubstrateOverlap2.eps", 
           width = 17.4/2.539998, height = 5.8/2.539998, pointsize = 10,
           family = 'TUM Neue Helvetica 55 Regular', encoding = 'Greek')

par(mar = c(3, 4, 5.5, 0.5), lwd = 0.01, xpd = TRUE)
bp <- barplot(sweep(rb, 2, colSums(rb), "/"), beside = FALSE,  ylab = "Proportion",
              border = "white", las = 2, lwd = 0.1, names.arg = rep("", ncol(rb)),
              col = .cs$colorScheme("c3"))

mtext(rb[1, ], line = 0.5, side = 3, at = bp, cex = 0.8, las = 2, col = .cs$colorScheme("c3")[1])
# mtext("Overlap", line = 0.5, side = 2, at = 1.2, las = 2, col = .cs$colorScheme("c3")[1])
mtext(rb[2, ], line = 1.5, side = 3, at = bp, cex = 0.8, las = 2, col = .cs$colorScheme("c3")[2])
mtext(rb[3, ], line = 3, side = 3, at = bp, cex = 0.8, las = 2, col = .cs$colorScheme("c3")[3])
mtext(colnames(rb), line = 0, side = 1, at = bp, cex = 0.8, las = 2)
legend(x = 0, y = 2.1, cex = 0.9, pt.cex = 1.5,
       legend = c("Recorded in database only", "Predicted by correlation only", "Overlap of the two"), 
       col = .cs$colorScheme("c3")[3:1], pch = 15, bty = "n", ncol=3)

# rownames(rb)
dev.off()

