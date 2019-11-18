# THIS IS WHAT IS USED FOR THE ATLANTIC SCORE!
setwd("/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen")
library(fastmatch)
library(GSVA)
library(matrixStats)

ev_prots_crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/ev_prots_crc65.rds')
pkinfam1 <- readRDS('/media/kusterlab/internal_projects/active/CRC65/data/24_CRC64_1and2and3_120ppm_1.4.1.2_peakproperties/analysis/pkinfam1.rda')
path <- readRDS(file = "Res/20180625_kinaseActivity/sub_Hallmarks_pathways.RDS") 
signor <- readRDS(file = "Res/20180625_kinaseActivity/sub_signor.RDS") 
integ <- readRDS(file = "Res/20180625_kinaseActivity/sub_integDB.RDS") 
nkin <- readRDS(file = "Res/20180625_kinaseActivity/sub_networkinSiteSub05.RDS") 
intnk <- readRDS(file = "Res/20180625_kinaseActivity/sub_integDB_networkin.RDS")

# newintegnames <- setdiff(rownames(data$CRC65$`Kinases (by Substrate phosphorylation)`), rownames(ev_prots_crc65))
newintegnames <- setdiff(names(integ), rownames(ev_prots_crc65))
newintegnames[nchar(newintegnames)<=2] <- paste0(newintegnames[nchar(newintegnames)<=2], '_EXT')

queryuniprot <- function(gene) {
  require(httr)
  require(jsonlite)
  require(xml2)
  gene <- gsub(' ', '%20', gene)
  cat(gene, '\n')
  requestURL <- paste0("https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&reviewed=true&gene=", gene, "&taxid=9606")
  r <- GET(requestURL, accept("application/xml"))
  stop_for_status(r)
  xml <- read_xml(content(r,as = "raw"))
  lst <- as_list(xml)
  nms <- unique(lapply(lst$uniprot, function(i) unlist(i$gene)))
  matches <- sapply(nms, function(i) any(i==gene))
  if (length(matches)==0) return('nomatch')
  res <- sapply(nms[matches], function(i) i[1])
  return(res)
}

potnames <- lapply(newintegnames, queryuniprot)
names(potnames) <- newintegnames
potnames[lapply(potnames, length)==0] <- 'nomatch'
mapped <- potnames[!unlist(lapply(potnames, function(i) all(i=='nomatch') | length(i)>1))]
unmapped <- potnames[unlist(lapply(potnames, function(i) all(i=='nomatch') | length(i)>1))]
unmapped[['AKT']] <- 'nomatch' # ambiguous
unmapped[['AMPK']] <- 'PRKAB1'
unmapped[['AURORA A']] <- 'AURKA'
unmapped[['AURORA B']] <- 'AURKB'
unmapped[['BCR/ABL']] <- 'ABL1'
unmapped[['BCR/ABL FUSION']] <- 'ABL1'
unmapped[['CAK']] <- 'nomatch' # ambiguous
unmapped[['CAM-KI_ALPHA']] <- 'CAMK1'
unmapped[['CAM-KII_ALPHA']] <- 'CAMK2A'
unmapped[['CAM-KIV']] <- 'CAMK4'
unmapped[['CAMK']] <- 'nomatch' # ambiguous
unmapped[['CDC42BP']] <- 'nomatch' # ambiguous
unmapped[['CK_EXT']] <- 'nomatch' # ambiguous
unmapped[['CK1']] <- 'CSNK1G2' # not exact match
unmapped[['CK1_ALPHA']] <- 'CSNK1A1'
unmapped[['CK1_DELTA']] <- 'CSNK1D'
unmapped[['CK1_EPSILON']] <- 'CSNK1E'
unmapped[['CK2']] <- 'nomatch' # ambiguous
unmapped[['CK2_ALPHA']] <- 'CSNK2A1'
unmapped[['CYCLINA2/CDK2']] <- 'CDK2'
unmapped[['CYCLINB/CDK1']] <- 'CDK1'
unmapped[['CYCLIND/CDK4']] <- 'CDK4'
unmapped[['CYCLINE/CDK2']] <- 'CDK2'
unmapped[['DNA-PK']] <- 'PRKDC'
unmapped[['EG3 KINASE']] <- 'MELK'
unmapped[['ERK1/2']] <- c('MAPK1', 'MAPK3')
unmapped[['GRK-2']] <- 'GRK2'
unmapped[['GRK-5']] <- 'GRK5'
unmapped[['GRK-6']] <- 'GRK6'
unmapped[['GSK-3_ALPHA']] <- 'GSK3A'
unmapped[['GSK-3_BETA']] <- 'GSK3B'
unmapped[['GSK3-ALPHA']] <- 'GSK3A'
unmapped[['GSK3-BETA']] <- 'GSK3B'
unmapped[['GSK3B/AXIN/APC']] <- 'GSK3B'
unmapped[['IKK_ALPHA']] <- 'CHUK'
unmapped[['IKK_BETA']] <- 'IKBKB'
unmapped[['MAK']] <- 'MAK' # from 2 databases - phophoNetwork and uniprot
unmapped[['MAPK1 AND MAPK3']] <- c('MAPK1', 'MAPK3')
unmapped[['MAST1']] <- 'MAST1'
unmapped[['MEK1/2']] <- c('MAP2K1', 'MAP2K2')
unmapped[['MRCKA']] <- 'CDC42BPA'
unmapped[['MTORC1']] <- 'MTOR'
unmapped[['MTORC2']] <- 'MTOR'
unmapped[['P38']] <- c('MAPK11', 'MAPK12', 'MAPK13', 'MAPK14')
unmapped[['P70S6K']] <- 'RPS6KB1'
unmapped[['P70S6KB']] <- 'RPS6KB2'
unmapped[['PAK']] <- c('PAK1', 'PAK2') # from 2 databases - signor and uniprot
unmapped[['PAK7/PAK5']] <- 'PAK5'
unmapped[['PDGFR_BETA']] <- 'PDGFRB'
unmapped[['PDK-1']] <- 'PDPK1'
unmapped[['PIM-1']] <- 'PIM1'
unmapped[['PKA']] <- c('PRKACA', 'PRKACB', 'PRKACG') # from 2 databases - signor and uniprot
unmapped[['PKA_ALPHA']] <- 'PRKACA'
unmapped[['PKB_BETA']] <- 'AKT2'
unmapped[['PKB/AKT1']] <- 'AKT1'
unmapped[['PKB/AKT2']] <- 'AKT2'
unmapped[['PKC']] <- 'PRKCZ'
unmapped[['PKC_ALPHA']] <- 'PRKCA'
unmapped[['PKC_BETA']] <- 'PRKCB'
unmapped[['PKC_DELTA']] <- 'PRKCD'
unmapped[['PKC_EPSILON']] <- 'PRKCE'
unmapped[['PKC_ETA']] <- 'PRKCH'
unmapped[['PKC_THETA']] <- 'PRKCQ'
unmapped[['PKC_ZETA']] <- 'PRKCZ'
unmapped[['PKC/PRKCA']] <- 'PRKCA'
unmapped[['PKC/PRKCB']] <- 'PRKCB'
unmapped[['PKC/PRKCD']] <- 'PRKCD'
unmapped[['PKC/PRKCE']] <- 'PRKCE'
unmapped[['PKC/PRKCG']] <- 'PRKCG'
unmapped[['PKC/PRKCH']] <- 'PRKCH'
unmapped[['PKC/PRKCI']] <- 'PRKCI'
unmapped[['PKC/PRKCQ']] <- 'PRKCQ'
unmapped[['PKC/PRKCZ']] <- 'PRKCZ'
unmapped[['PKD/PRKD1']] <- 'PRKD1'
unmapped[['PKD/PRKD2']] <- 'PRKD2'
unmapped[['PKD1']] <- 'PRKD1'
unmapped[['PKG/PRKG1']] <- 'PRKG1'
unmapped[['PKG1/CGK-I']] <- 'PRKG1'
unmapped[['PKG2/CGK-II']] <- 'PRKG2'
unmapped[['RAF']] <- 'RAF1'
unmapped[['RAGE']] <- 'MOK'
unmapped[['RSK-1']] <- 'RPS6KA1'
unmapped[['RSK-2']] <- 'RPS6KA2'
unmapped[['STK3/MST2']] <- 'STK3'
unmapped[['STK4/MST1']] <- 'STK4'
unmapped[['ZIPK/DAPK3']] <- 'DAPK3'
mapped2 <- unmapped[!unlist(lapply(unmapped, function(i) all(i=='nomatch') | length(i)>1))]
unmapped <- unmapped[unlist(lapply(unmapped, function(i) all(i=='nomatch') | length(i)>1))]
for (i in seq_along(unmapped)) {
  unmapped[[i]] <- names(unmapped)[i]
}
for (i in seq_along(mapped)) {
  names(mapped[[i]]) <- NULL
}
for (i in seq_along(mapped2)) {
  names(mapped2[[i]]) <- NULL
}

mapping <- unlist(c(mapped, mapped2, unmapped))


for (i in newintegnames) {
  idx <- match(i, names(integ))
  names(integ)[idx] <- mapping[i]
}
newinteg <- list()
for (i in unique(names(integ))) {
  newinteg[[i]] <- unlist(integ[names(integ)%in%i])
  names(newinteg[[i]]) <- NULL
}
integ <- newinteg
integ <- integ[names(integ)%in%pkinfam1$Name]

require(data.table)
mappingdt <- data.table(`Kinase accession`=names(mapping), `Kinase gene name`=unname(mapping))
mappingdt <- mappingdt[`Kinase gene name`%in%pkinfam1$Name,]
mappingdt <- rbindlist(list(mappingdt, data.table(`Kinase accession`=setdiff(names(integ),mappingdt$`Kinase accession`),
                                                  `Kinase gene name`=setdiff(names(integ),mappingdt$`Kinase accession`))))
fwrite(mappingdt, file = '/media/kusterlab/internal_projects/active/NCI60_Phospho/Manuscript/current/supplementTables/new/ksmapping.txt', sep = '\t', col.names = T, row.names = F)

# 
# ##### GSVA #####
tryp <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")
emat <- tryp$fun_prepElnet(tryp, panel = "NCI60")
emat <- emat$imputed

score_nkin <- gsva(emat, nkin, method = "zscore")
score_path <- gsva(emat, path, method = "zscore")
score_integ  <- gsva(emat, integ, method = "zscore")
score_signor <- gsva(emat, signor, method = "zscore")
score_intnk <- gsva(emat, intnk, method = "zscore")

scorelist <- list(nkin = score_nkin, 
                  path = score_path, 
                  integ = score_integ, 
                  signor = score_signor,
                  integ_nkin = score_intnk
                  )

saveRDS(scorelist, file = "Res/20180625_kinaseActivity/scorelist_nci60.RDS")

##3
emat <- tryp$fun_prepElnet(tryp, panel = "CRC65", celllines = setdiff(colnames(tryp$crc65_expr), "CoCM-1_CRC65"))
emat <- emat$imputed

score_nkin <- gsva(emat, nkin, method = "zscore")
score_path <- gsva(emat, path, method = "zscore")
score_integ  <- gsva(emat, integ, method = "zscore")
score_signor <- gsva(emat, signor, method = "zscore")
score_intnk <- gsva(emat, intnk, method = "zscore")

scorelist <- list(nkin = score_nkin, 
                  path = score_path, 
                  integ = score_integ, 
                  signor = score_signor,
                  integ_nkin = score_intnk)

saveRDS(scorelist, file = "Res/20180625_kinaseActivity/scorelist_crc65.RDS")

