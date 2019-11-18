source('Shiny/prepShiny.R')

##
sum1 <- res.elnet$EFFECTSIZE * res.elnet$FREQ
sum1split <- split(sum1, paste(res.elnet$DRUG, res.elnet$INPUTDATA))

top100 <- lapply(sum1split, function(x) {
  rk <- rank(-abs(x))
  rk < 100
})
top100 <- unlist(top100)
freq50 <- res.elnet$FREQ > 50

i <- top100 & freq50
table(i)

sub.elnet <- res.elnet[i, ]

ll <- apply(sub.elnet, 1, function(x) {
  getAnnot(gene = x["ID"], inputdata = x["INPUTDATA"], col = a.cols)
})
llmat <- do.call(rbind, ll)
mat <- cbind(sub.elnet, llmat)

write.table(mat, file = "Res/20170112_elnetResSummary/elnetSummary.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")





