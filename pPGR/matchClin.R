library(openxlsx)

cli <- read.xlsx("Dat/PGR/Silver Kollektiv Follow up_survivalstatus.xlsx")
cli2 <- read.xlsx("../Progesteronrezeptor/Kopie von Silver Kollektiv Follow up.xlsx")
identical(cli$TMA.number, cli2$`Int..#`[1:363])
cli$Patho <- cli2$`Patho#`[1:363]


pgr <- readRDS("Res/20170707_pGR/pgrstain.RDS")
pgr$stain.intensity.phos <- as.numeric(pgr$stain.intensity.phos)
head(pgr)


blk <- as.character(pgr$HNR.block)
blk <- strsplit(blk, " ")
blk <- sapply(blk, "[", 1)
pgr$cid <- blk
pgr$h.protein <- pgr$stain.percent.protein * pgr$stain.intensity.protein
pgr$h.phos <- pgr$stain.percent.phos * pgr$stain.intensity.phos

path <- as.character(cli$Patho)
path <- sapply(strsplit(path, " "), "[", 1)
cli$cid <- path


#### 
x <- tapply(pgr$h.protein, pgr$cid, mean, na.rm = TRUE)
cli$h.score.mean.protein <- x[cli$cid]

x <- tapply(pgr$h.protein, pgr$cid, max, na.rm = TRUE)
cli$h.score.max.protein <- x[cli$cid]

x <- tapply(pgr$h.phos, pgr$cid, mean, na.rm = TRUE)
cli$h.score.mean.phos <- x[cli$cid]

x <- tapply(pgr$h.phos, pgr$cid, max, na.rm = TRUE)
cli$h.score.max.phos <- x[cli$cid]

#
x <- tapply(pgr$stain.percent.protein, pgr$cid, mean, na.rm = TRUE)
cli$sperc.mean.protein <- x[cli$cid]

x <- tapply(pgr$stain.percent.protein, pgr$cid, max, na.rm = TRUE)
cli$sperc.max.protein <- x[cli$cid]

x <- tapply(pgr$stain.percent.phos, pgr$cid, mean, na.rm = TRUE)
cli$sperc.mean.phos <- x[cli$cid]

x <- tapply(pgr$stain.percent.phos, pgr$cid, max, na.rm = TRUE)
cli$sperc.max.phos <- x[cli$cid]


# sintensi.mean.protein
x <- tapply(pgr$stain.intensity.protein, pgr$cid, mean, na.rm = TRUE)
cli$sintens.mean.protein <- x[cli$cid]

x <- tapply(pgr$stain.intensity.protein, pgr$cid, max, na.rm = TRUE)
cli$sintens.max.protein <- x[cli$cid]

x <- tapply(pgr$stain.intensity.phos, pgr$cid, mean, na.rm = TRUE)
cli$sintens.mean.phos <- x[cli$cid]

x <- tapply(pgr$stain.intensity.phos, pgr$cid, max, na.rm = TRUE)
cli$sintens.max.phos <- x[cli$cid]



for (i in names(cli)) 
  cli[[i]][is.infinite(cli[[i]])] <- NA

#
saveRDS(cli, file = "Res/20170707_pGR/pgrstainReady.RDS")







