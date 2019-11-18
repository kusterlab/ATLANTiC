if (!exists("matchobject")) {
  matchobject <- eval(as.symbol(as.character(readline(prompt = 'Match annotation to which vector of cell line names (e.g. tt after tt <- colnames(log10lfqkmri)) or cells?' ))))
}

# Bracht & Bodmer
FUannotation <- read.table("/media/kusterlab/internal_projects/active/CRC65/data/24_CRC64_1and2and3_120ppm_1.4.1.2_peakproperties/analysis/5-Fluorouracil response in a large panel of co.csv", header=T, sep=",", stringsAsFactors=F)
FUannotation <- FUannotation[c(1:(dim(FUannotation)[1]-3)),]
FUannotation$Cell.Line <- sapply(FUannotation$Cell.Line, function(dat) gsub(pattern="[\t]*", "", x=dat))
FUannotation[which(FUannotation$Cell.Line=="CACO2"),1] <- "CaCo-2"
FUannotation[which(FUannotation$Cell.Line=="CCK81"),1] <- "CCK-81"
FUannotation[which(FUannotation$Cell.Line=="COLO320DM"),1] <- "Colo 320DM"
FUannotation[which(FUannotation$Cell.Line=="C125PM"),1] <- "C125-PM"
FUannotation[which(FUannotation$Cell.Line=="CAR1"),1] <- "CaR-1"
FUannotation[which(FUannotation$Cell.Line=="CoCM1"),1] <- "CoCM-1"
FUannotation[which(FUannotation$Cell.Line=="COLO678"),1] <- "Colo-678"
FUannotation[which(FUannotation$Cell.Line=="COLO741"),1] <- "Colo 741"
FUannotation[which(FUannotation$Cell.Line=="DLD1"),1] <- "DLD-1"
FUannotation[which(FUannotation$Cell.Line=="HDC111"),1] <- "HDC-111"
FUannotation[which(FUannotation$Cell.Line=="HDC135"),1] <- "HDC-135"
FUannotation[which(FUannotation$Cell.Line=="HDC142"),1] <- "HDC-142"
FUannotation[which(FUannotation$Cell.Line=="HDC54"),1] <- "HDC-54"
FUannotation[which(FUannotation$Cell.Line=="HDC57"),1] <- "HDC-57"
FUannotation[which(FUannotation$Cell.Line=="HDC73"),1] <- "HDC-73"
FUannotation[which(FUannotation$Cell.Line=="HDC8"),1] <- "HDC-8"
FUannotation[which(FUannotation$Cell.Line=="HDC82"),1] <- "HDC-82"
FUannotation[which(FUannotation$Cell.Line=="HDC9"),1] <- "HDC-9"
FUannotation[which(FUannotation$Cell.Line=="RCM1"),1] <- "RCM-1"
FUannotation[which(FUannotation$Cell.Line=="HCA46"),1] <- "HCA-46"
FUannotation[which(FUannotation$Cell.Line=="HCA7"),1] <- "HCA-7"
FUannotation[which(FUannotation$Cell.Line=="HCT116"),1] <- "HCT 116"
FUannotation[which(FUannotation$Cell.Line=="HRA19"),1] <- "HRA-19"
FUannotation[which(FUannotation$Cell.Line=="HT29"),1] <- "HT-29"
FUannotation[which(FUannotation$Cell.Line=="LOVO"),1] <- "LoVo"
FUannotation[which(FUannotation$Cell.Line=="LS174T"),1] <- "LS 174T"
FUannotation[which(FUannotation$Cell.Line=="LS180"),1] <- "LS 180"
FUannotation[which(FUannotation$Cell.Line=="NCIH548"),1] <- "NCI-H548"
FUannotation[which(FUannotation$Cell.Line=="NCIH716"),1] <- "NCI-H716"
FUannotation[which(FUannotation$Cell.Line=="NCIH747"),1] <- "NCI-H747"
FUannotation[which(FUannotation$Cell.Line=="OXCO1"),1] <- "OXCO-1"
FUannotation[which(FUannotation$Cell.Line=="OXCO3"),1] <- "OXCO-3"
FUannotation[which(FUannotation$Cell.Line=="SKCO1"),1] <- "SK-CO-1"
FUannotation[which(FUannotation$Cell.Line=="SNUC2B"),1] <- "SNU-C2B"
FUannotation[which(FUannotation$Cell.Line=="VACO4A"),1] <- "VACO 4A"
rownames(FUannotation) <- FUannotation$Cell.Line
# FUannotation <- FUannotation[intersect(rownames(cormatttest2), rownames(FUannotation)),]
# FUannotation <- merge(FUannotation, data.frame(Cell.Line=colnames(matchobject)), by="Cell.Line", all=T)
FUannotation <- FUannotation[intersect(matchobject, rownames(FUannotation)),]
FUannotation <- merge(FUannotation, data.frame(Cell.Line=matchobject), by="Cell.Line", all=T)
rownames(FUannotation) <- FUannotation$Cell.Line
# FUannotation <- FUannotation[match(colnames(cormatttest2), rownames(FUannotation)),]
FUannotation <- FUannotation[match(matchobject, rownames(FUannotation)),]
FUannotation$Cell.Line <- NULL
FUannotation$Passage.by <- NULL
FUannotation$Media <- NULL
FUannotation$Source <- NULL
FUannotation$Reference..RER.status. <- NULL
FUannotation$Lag.time..days. <- as.double(FUannotation$Lag.time..days.)
FUannotation$Doubling.time..hours. <- as.double(FUannotation$Doubling.time..hours.)
colnames(FUannotation)[grep("Lag.time..days.", colnames(FUannotation))] <- "Lag time (days)"
colnames(FUannotation)[grep("Doubling.time..hours.", colnames(FUannotation))] <- "Doubling time (hours)"
colnames(FUannotation) <- sapply(colnames(FUannotation), function(dat) gsub(pattern="\\.", replacement=" ", dat))