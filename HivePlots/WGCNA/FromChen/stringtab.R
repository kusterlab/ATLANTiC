library(data.table)
library(fastmatch)
library(stringr)
library(reshape2)

################################### processing bioplex ##################################
bioplex <- fread("dat/BIOPLEX/bioplex20_hsa", integer64 = "double", data.table = FALSE)
head(bioplex)
hist(bioplex$`p(Interaction)`)

bp <- data.frame(p1 = pmin(bioplex$SymbolA, bioplex$SymbolB), 
                 p2 = pmax(bioplex$SymbolA, bioplex$SymbolB), 
                 stringsAsFactors = FALSE)
bp <- bp[bp$p1 != bp$p2, ]
bp0 <- paste(bp$p1, bp$p2)



################################### processing string ##################################
str <- fread("Dat/STRING/9606.protein.links.detailed.v10.5.txt", integer64 = "double", data.table = FALSE)
hist(str$combined_score)
str <- str[str$combined_score > 750, ]
str$protein1 <- gsub("9606.", "", str$protein1)
str$protein2 <- gsub("9606.", "", str$protein2)
str <- str[order(str$combined_score, decreasing = TRUE), ]
idmap <- fread("Dat/Uniprot/HUMAN_9606_idmapping.dat", data.table = FALSE, header = FALSE)
str$uniprotA <- idmap[fmatch(str$protein1, idmap$V3), 1]
str$uniprotB <- idmap[fmatch(str$protein2, idmap$V3), 1]
str$uniprotA <- str_split_fixed(string = str$uniprotA, pattern = "-", n = 2)[, 1]
str$uniprotB <- str_split_fixed(string = str$uniprotB, pattern = "-", n = 2)[, 1]
pt <- paste(idmap$V1, idmap$V2, sep = "_")
str$symbolA <- idmap[fmatch(paste(str$uniprotA, "Gene_Name", sep = "_"), pt), 3]
str$symbolB <- idmap[fmatch(paste(str$uniprotB, "Gene_Name", sep = "_"), pt), 3]
str <- str[! (is.na(str$symbolA) | is.na(str$symbolB)), ]
str0 <- data.frame(p1 = pmin(str$symbolA, str$symbolB), 
                   p2 = pmax(str$symbolA, str$symbolB), 
                   stringsAsFactors = FALSE)
str0 <- str0[str0$p1 != str0$p2, ]
tstr0 <- paste(str0$p1, str0$p2)
tstr0 <- unique(tstr0)
gc()


################################### processing site ppi ##################################
ppicrc <- readRDS("Res/20180314_wgcnaSupTables/crc65_cormat.RDS")
ppicrc$Var1 <- str_split_fixed(as.character(ppicrc$Var1), "_", 2)[, 1]
ppicrc$Var2 <- str_split_fixed(as.character(ppicrc$Var2), "_", 2)[, 1]
crc <- data.frame(p1 = pmin(ppicrc$Var1, ppicrc$Var2), 
                  p2 = pmax(ppicrc$Var1, ppicrc$Var2), 
                  cor = ppicrc$value, 
                  stringsAsFactors = FALSE)
crc$label <- paste(crc$p1, crc$p2)

ppinci <- readRDS("Res/20180314_wgcnaSupTables/nci60_cormat.RDS")
ppinci$Var1 <-  str_split_fixed(as.character(ppinci$Var1), "_", 2)[, 1]
ppinci$Var2 <- str_split_fixed(as.character(ppinci$Var2), "_", 2)[, 1]
nci <- data.frame(p1 = pmin(ppinci$Var1, ppinci$Var2), 
                  p2 = pmax(ppinci$Var1, ppinci$Var2), 
                  cor = ppinci$value, 
                  stringsAsFactors = FALSE)
nci$label <- paste(nci$p1, nci$p2)







library(VennDiagram)
vv <- sapply(c(0.5,0.6, 0.7, 0.8, 0.9), function(cut) {
  print(cut)
  q1 <- nci[ abs(nci$cor) > cut, ]
  q2 <- crc[ abs(crc$cor) > cut, ]
  
  s1 <- unique(q1$label[q1$p1 != q1$p2])
  s2 <- unique(q2$label[q2$p1 != q2$p2])
  
  ll <- list(NCI60 = s1, CRC65 = s2, STRING = tstr0, bioplex = bp0)
  
  file <- file.path("Res/20180314_wgcnaSupTables/", paste0("corcut", cut, ".png"))
  venn.diagram(ll, filename = file)
})

