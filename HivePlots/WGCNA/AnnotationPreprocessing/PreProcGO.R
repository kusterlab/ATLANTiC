require(org.At.tair.db)
require(org.Hs.eg.db)
require(org.Mm.eg.db)
require(GO.db)
require(reactome.db)
require(parallel)
require(data.table)
require(reshape2)
require(RCurl)

source('src/download_annotation.R')

# ---------------------------------------- #
# Arabidopsis thaliana
# ---------------------------------------- #

ath_idmap <- data.table(select(x = org.At.tair.db, keys = keys(org.At.tair.db, keytype = 'TAIR'), columns=c('ENTREZID', 'SYMBOL'), keytype = 'TAIR'))
invisible(lapply(colnames(ath_idmap), function(i) ath_idmap[eval(as.name(i))=='',eval(i):=NA]))

# ath_symbol <- ath_idmap[,list(SYMBOL=paste0(sort(unique(SYMBOL)), collapse = ';')),by=TAIR]
# setkey(ath_symbol, TAIR)
# ath_entrez <- ath_idmap[,list(ENTREZID=paste0(sort(unique(ENTREZID)), collapse = ';')),by=TAIR]
# setkey(ath_entrez, TAIR)
# ath_entrez[ENTREZID=='',ENTREZID:=NA]

ath_go <- fread(
  'dat/annotations/ath_go',
  stringsAsFactors = F,
  integer64 = 'double',
  col.names = c(
    'DB',
    'DB_Object_ID',
    'DB_Object_Symbol',
    'Qualifier (optional)',
    'GO ID',
    'DB:Reference(|DB:Reference)',
    'Evidence',
    'With (or) From (optional)',
    'Aspect',
    'DB_Object_Name(|Name) (optional)',
    'DB_Object_Synonym(|Synonym) (optional)',
    'DB_Object_Type',
    'taxon(|taxon)',
    'Date',
    'Assigned_by',
    'V1',
    'DB:DB_Object_ID')
)
ath_go[,V1:=NULL]
ath_go[,ID:=`GO ID`]
ath_go[,TAIR:=`DB_Object_Name(|Name) (optional)`]
ath_go[,SYMBOL:=DB_Object_Symbol]
invisible(ath_go[Aspect=='P',Aspect:='BP'])
invisible(ath_go[Aspect=='F',Aspect:='MF'])
invisible(ath_go[Aspect=='C',Aspect:='CC'])
ath_go[,ONTOLOGY:=paste0('GO', '_', Aspect)]
ath_go <- ath_go[,list(TAIR, SYMBOL, ID, ONTOLOGY, present = T)]
invisible(lapply(colnames(ath_go), function(i) ath_go[eval(as.name(i))=='',eval(i):=NA]))

setkey(ath_idmap, TAIR)
setkey(ath_go, TAIR)
ath_go <- merge(ath_go, ath_idmap, by='TAIR', allow.cartesian = T, all = T)
ath_go[is.na(SYMBOL.x),SYMBOL.x:=SYMBOL.y]
ath_go[,SYMBOL.y:=NULL]
setnames(ath_go, 'SYMBOL.x', 'SYMBOL')

ath_cor <- fread('/media/kusterlab/internal_projects/active/Arabidopsis_proteome/Arabdiopsis thaliana/SEC/ortholog complex list_McBride2017_JM.txt')
setnames(ath_cor, c('Complex name', 'Complex number', 'ATG Arabidopsis'), c('Name', 'ID', 'TAIR'))
ath_cor <- ath_cor[,list(ID,
                         Name,
                         TAIR)]
setkey(ath_idmap, TAIR)
setkey(ath_cor, TAIR)
ath_cor <- merge(ath_cor, ath_idmap, by='TAIR', allow.cartesian = T, all.x = T)

ath_cor[,ONTOLOGY:='CORUM']
ath_cor <- ath_cor[!is.na(ID)]
ath_cor[,present:=T]
ath_cor[,ID:=paste0('CORUM_ath_', ID)]
ath_cor_names <- ath_cor[,list(ID, Name)]
setkey(ath_cor_names)
ath_cor_names <- unique(ath_cor_names)
ath_cor <- ath_cor[,list(TAIR, SYMBOL, ENTREZID, ID, ONTOLOGY, present)]
setkey(ath_cor)
ath_cor <- unique(ath_cor)

ath_kegg <- fread('dat/annotations/keggedges_ath', header = F)
setnames(ath_kegg, c('ID','TAIR'))
ath_kegg[,ID:=gsub('^path\\:', '', ID)]
ath_kegg[,TAIR:=gsub('^ath\\:', '', TAIR)]
setcolorder(ath_kegg, c('TAIR', 'ID'))
invisible(lapply(colnames(ath_kegg), function(i) ath_kegg[eval(as.name(i))=='',eval(i):=NA]))

setkey(ath_idmap, TAIR)
setkey(ath_kegg, TAIR)
ath_kegg <- merge(ath_kegg, ath_idmap, by='TAIR', allow.cartesian = T, all = T)
# ath_kegg[is.na(SYMBOL.x),SYMBOL.x:=SYMBOL.y]
# ath_kegg[,SYMBOL.y:=NULL]
# setnames(ath_go, 'SYMBOL.x', 'SYMBOL')
ath_kegg[,ONTOLOGY:='KEGG']
ath_kegg[,present:=T]

ath_rct <- fread(
  'dat/annotations/reactomepathways_ath',
  stringsAsFactors = F,
  integer64 = 'double',
  col.names = c(
    'Source_database_identifier',
    'Reactome_Stable_identifier',
    'URL',
    'Event_(Pathway_or_Reaction)_Name',
    'Evidence_Code',
    'Species'
  )
)

ath_rct <- ath_rct[Species=='Arabidopsis thaliana',]

setnames(ath_rct, c('Source_database_identifier', 'Reactome_Stable_identifier','Event_(Pathway_or_Reaction)_Name'), c('TAIR','ID','Name'))
ath_rct <- ath_rct[,list(TAIR, ID, Name)]
invisible(lapply(colnames(ath_rct), function(i) ath_rct[eval(as.name(i))=='',eval(i):=NA]))

setkey(ath_idmap, TAIR)
setkey(ath_rct, TAIR)
ath_rct <- merge(ath_rct, ath_idmap, by='TAIR', allow.cartesian = T, all.x = T)
ath_rct[,ONTOLOGY:='REACTOME']
ath_rct[,present:=T]
ath_rct_names <- ath_rct[,list(ID, Name)]
setkey(ath_rct_names)
ath_rct_names <- unique(ath_rct_names)
ath_rct <- ath_rct[,list(TAIR, SYMBOL, ENTREZID, ID, ONTOLOGY, present)]
setkey(ath_rct)
ath_rct <- unique(ath_rct)

ath <- rbindlist(list(ath_go, ath_kegg, ath_rct, ath_cor), use.names = T)[!is.na(present),]

setkey(ath)
ath <- unique(ath)
#
# setkey(ath, ID)
# ath_offspring <- rbindlist(
#   list(
#     ath,
#     lapply(ath[grepl('^GO_', ONTOLOGY),unique(ONTOLOGY)][1], function(i)
#       rbindlist(
#         mclapply(ath[ONTOLOGY==i,unique(ID)][1:10], function(j)
#           rbindlist(
#             lapply(
#                , function(k)
#               {
#                 tmp <- ath[.(j),]
#                 tmp[,ID:=k]
#                 return(tmp)
#               }
#             )
#           ),
#           mc.cores = 20
#         )
#       )
#     )
#   )
# )


# # PREVIOUS VERSION!
# # sapply(columns(org.At.tair.db), function(i) head(keys(org.At.tair.db, keytype = i)))
# ath <- data.table(select(x = org.At.tair.db, keys = keys(org.At.tair.db, keytype = 'TAIR'), columns=c('ENTREZID', 'SYMBOL', 'PATH', 'GO'), keytype = 'TAIR'))
# ath[,EVIDENCE:=NULL]
# ath <- ath[!is.na(PATH) | !is.na(GO) | !is.na(ONTOLOGY),]
# setkey(ath)
# ath <- unique(ath)
# setnames(ath, 'GO', 'ID')
# ath <- rbindlist(list(ath[!is.na(PATH),list(TAIR,ENTREZID,SYMBOL,ID=PATH,ONTOLOGY='KEGG')],ath[!is.na(ID),list(TAIR,ENTREZID,SYMBOL,ID,ONTOLOGY=paste0('GO_', ONTOLOGY))]))
# setkey(ath)
# ath <- unique(ath)
# ath[ONTOLOGY=='KEGG',ID:=paste0('ath',ID)]
# ath[,present:=T]
# # ath_func <- castfun('TAIR',ath)

# ---------------------------------------- #
# Homo sapiens
# ---------------------------------------- #

hsa_idmap <- data.table(select(x = org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = 'UNIPROT'), columns=c('ENTREZID', 'SYMBOL'), keytype = 'UNIPROT'))
invisible(lapply(colnames(hsa_idmap), function(i) hsa_idmap[eval(as.name(i))=='',eval(i):=NA]))

# hsa_symbol <- hsa_idmap[,list(SYMBOL=paste0(sort(unique(SYMBOL)), collapse = ';')),by=TAIR]
# setkey(hsa_symbol, TAIR)
# hsa_entrez <- hsa_idmap[,list(ENTREZID=paste0(sort(unique(ENTREZID)), collapse = ';')),by=TAIR]
# setkey(hsa_entrez, TAIR)
# hsa_entrez[ENTREZID=='',ENTREZID:=NA]

hsa_go <- fread(
  'dat/annotations/hsa_go',
  stringsAsFactors = F,
  integer64 = 'double',
  col.names = c(
    'DB',
    'DB_Object_ID',
    'DB_Object_Symbol',
    'Qualifier (optional)',
    'GO ID',
    'DB:Reference(|DB:Reference)',
    'Evidence',
    'With (or) From (optional)',
    'Aspect',
    'DB_Object_Name(|Name) (optional)',
    'DB_Object_Synonym(|Synonym) (optional)',
    'DB_Object_Type',
    'taxon(|taxon)',
    'Date',
    'Assigned_by',
    'Annotation_Extension (optional)',
    'Gene_Product_Form_ID (optional)')
)
# hsa_go[,V1:=NULL]
hsa_go[,ID:=`GO ID`]
hsa_go[,UNIPROT:=`DB_Object_ID`]
hsa_go[,SYMBOL:=DB_Object_Symbol]
invisible(hsa_go[Aspect=='P',Aspect:='BP'])
invisible(hsa_go[Aspect=='F',Aspect:='MF'])
invisible(hsa_go[Aspect=='C',Aspect:='CC'])
hsa_go[,ONTOLOGY:=paste0('GO', '_', Aspect)]
hsa_go <- hsa_go[,list(UNIPROT, SYMBOL, ID, ONTOLOGY, present = T)]
invisible(lapply(colnames(hsa_go), function(i) hsa_go[eval(as.name(i))=='',eval(i):=NA]))

setkey(hsa_idmap, UNIPROT)
setkey(hsa_go, UNIPROT)
hsa_go <- merge(hsa_go, hsa_idmap, by='UNIPROT', allow.cartesian = T, all = T)
hsa_go[is.na(SYMBOL.x),SYMBOL.x:=SYMBOL.y]
hsa_go[,SYMBOL.y:=NULL]
setnames(hsa_go, 'SYMBOL.x', 'SYMBOL')
hsa_go <- hsa_go[!is.na(ID),]

hsa_cor <- fread('dat/annotations/coreComplexes.txt', stringsAsFactors = F, integer64 = 'double')
hsa_cor <- hsa_cor[Organism=='Human',]
hsa_cor <- hsa_cor[,list(ComplexID,
                         ComplexName,
                         `subunits(UniProt IDs)`,
                         `subunits(Entrez IDs)`,
                         `subunits(Gene name)`)]
hsa_cor <- hsa_cor[,list(Name=ComplexName,
                         UNIPROT=unique(unlist(strsplit(`subunits(UniProt IDs)`, ';')))),by=ComplexID]
setnames(hsa_cor, 'ComplexID', 'ID')
# UNIPROT=unique(unlist(strsplit(`subunits(UniProt IDs)`, ';'))),
# SYMBOL=unique(unlist(strsplit(`subunits(Gene name)`, ';'))),
setkey(hsa_idmap, UNIPROT)
setkey(hsa_cor, UNIPROT)
hsa_cor <- merge(hsa_cor, hsa_idmap, by='UNIPROT', allow.cartesian = T, all.x = T)
# hsa_kegg[is.na(SYMBOL.x),SYMBOL.x:=SYMBOL.y]
# hsa_kegg[,SYMBOL.y:=NULL]
# setnames(hsa_go, 'SYMBOL.x', 'SYMBOL')
hsa_cor[,ONTOLOGY:='CORUM']
hsa_cor <- hsa_cor[!is.na(ID)]
hsa_cor[,present:=T]
hsa_cor[,ID:=paste0('CORUM_hsa_', ID)]
hsa_cor_names <- hsa_cor[,list(ID, Name)]
setkey(hsa_cor_names)
hsa_cor_names <- unique(hsa_cor_names)
hsa_cor <- hsa_cor[,list(UNIPROT, SYMBOL, ENTREZID, ID, ONTOLOGY, present)]
setkey(hsa_cor)
hsa_cor <- unique(hsa_cor)

hsa_kegg <- fread('dat/annotations/keggedges_hsa', header = F)
setnames(hsa_kegg, c('ID','ENTREZID'))
hsa_kegg[,ID:=gsub('^path\\:', '', ID)]
hsa_kegg[,ENTREZID:=gsub('^hsa\\:', '', ENTREZID)]
setcolorder(hsa_kegg, c('ENTREZID', 'ID'))
invisible(lapply(colnames(hsa_kegg), function(i) hsa_kegg[eval(as.name(i))=='',eval(i):=NA]))

setkey(hsa_idmap, ENTREZID)
setkey(hsa_kegg, ENTREZID)
hsa_kegg <- merge(hsa_kegg, hsa_idmap, by='ENTREZID', allow.cartesian = T, all = T)
# hsa_kegg[is.na(SYMBOL.x),SYMBOL.x:=SYMBOL.y]
# hsa_kegg[,SYMBOL.y:=NULL]
# setnames(hsa_go, 'SYMBOL.x', 'SYMBOL')
hsa_kegg[,ONTOLOGY:='KEGG']
hsa_kegg <- hsa_kegg[!is.na(ID)]
hsa_kegg[,present:=T]


hsa_rct <- fread(
  'dat/annotations/reactomepathways_hsa',
  stringsAsFactors = F,
  integer64 = 'double',
  col.names = c(
    'Source_database_identifier',
    'Reactome_Stable_identifier',
    'URL',
    'Event_(Phsaway_or_Reaction)_Name',
    'Evidence_Code',
    'Species'
  )
)

hsa_rct <- hsa_rct[Species=='Homo sapiens',]

setnames(hsa_rct, c('Source_database_identifier', 'Reactome_Stable_identifier','Event_(Phsaway_or_Reaction)_Name'), c('UNIPROT','ID','Name'))
hsa_rct <- hsa_rct[,list(UNIPROT, ID, Name)]
invisible(lapply(colnames(hsa_rct), function(i) hsa_rct[eval(as.name(i))=='',eval(i):=NA]))

setkey(hsa_idmap, UNIPROT)
setkey(hsa_rct, UNIPROT)
hsa_rct <- merge(hsa_rct, hsa_idmap, by='UNIPROT', allow.cartesian = T, all.x = T)
hsa_rct[,ONTOLOGY:='REACTOME']
hsa_rct[,present:=T]
hsa_rct_names <- hsa_rct[,list(ID, Name)]
setkey(hsa_rct_names)
hsa_rct_names <- unique(hsa_rct_names)
hsa_rct <- hsa_rct[,list(UNIPROT, SYMBOL, ENTREZID, ID, ONTOLOGY, present)]
setkey(hsa_rct)
hsa_rct <- unique(hsa_rct)


hsa <- rbindlist(list(hsa_go, hsa_kegg, hsa_rct, hsa_cor), use.names = T)[!is.na(present),]

setkey(hsa)
hsa <- unique(hsa)



# # PREVIOUS VERSION
# # sapply(columns(org.Hs.eg.db), function(i) head(keys(org.Hs.eg.db, keytype = i)))
# hsa <- data.table(select(x = org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = 'UNIPROT'), columns=c('ENTREZID', 'SYMBOL', 'PATH', 'GO'), keytype = 'UNIPROT'))
# hsa[,EVIDENCE:=NULL]
# hsa <- hsa[!is.na(PATH) | !is.na(GO) | !is.na(ONTOLOGY),]
# setkey(hsa)
# hsa <- unique(hsa)
# setnames(hsa, 'GO', 'ID')
# hsa <- rbindlist(list(hsa[!is.na(PATH),list(UNIPROT,ENTREZID,SYMBOL,ID=PATH,ONTOLOGY='KEGG')],hsa[!is.na(ID),list(UNIPROT,ENTREZID,SYMBOL,ID,ONTOLOGY=paste0('GO_', ONTOLOGY))]))
# setkey(hsa)
# hsa <- unique(hsa)
# hsa[ONTOLOGY=='KEGG',ID:=paste0('hsa',ID)]
# hsa[,present:=T]
# # hsa_func <- castfun('UNIPROT',hsa)
# # sapply(columns(reactome.db), function(i) head(keys(reactome.db, keytype = i)))
# rct <- data.table(select(x = reactome.db, keys = keys(reactome.db, keytype = 'ENTREZID'), columns=c('ENTREZID', 'PATHID', 'PATHNAME'), keytype = 'ENTREZID'))
# rct[,org:=unlist(lapply(strsplit(PATHNAME, ': '), function(i) i[[1]]))]
# rct[,PATHNAME:=unlist(lapply(strsplit(PATHNAME, ': '), function(i) i[[2]]))]
# rct <- rct[org%in%c('Homo sapiens'),]
# setnames(rct, colnames(rct), c('ENTREZID', 'ID', 'Name', 'org'))
# rct[,org:=NULL]
# rct[,ONTOLOGY:='Reactome']
# 
# mrg <- merge(rct, hsa[,list(ENTREZID,UNIPROT,SYMBOL)], by='ENTREZID', all.x = T, allow.cartesian = T)
# setkey(mrg)
# mrg <- unique(mrg)
# mrg[,present:=T]
# 
# hsa <- rbindlist(list(hsa, mrg[,list(UNIPROT, ENTREZID, SYMBOL, ID, ONTOLOGY, present)]))


# ---------------------------------------- #
# Mus musculus
# ---------------------------------------- #

mmu_idmap <- data.table(select(x = org.Mm.eg.db, keys = keys(org.Mm.eg.db, keytype = 'UNIPROT'), columns=c('ENTREZID', 'SYMBOL'), keytype = 'UNIPROT'))
invisible(lapply(colnames(mmu_idmap), function(i) mmu_idmap[eval(as.name(i))=='',eval(i):=NA]))

# mmu_symbol <- mmu_idmap[,list(SYMBOL=paste0(sort(unique(SYMBOL)), collapse = ';')),by=TAIR]
# setkey(mmu_symbol, TAIR)
# mmu_entrez <- mmu_idmap[,list(ENTREZID=paste0(sort(unique(ENTREZID)), collapse = ';')),by=TAIR]
# setkey(mmu_entrez, TAIR)
# mmu_entrez[ENTREZID=='',ENTREZID:=NA]

mmu_go <- fread(
  'dat/annotations/mmu_go',
  skip = 46,
  stringsAsFactors = F,
  integer64 = 'double',
  col.names = c(
    'DB',
    'DB_Object_ID',
    'DB_Object_Symbol',
    'Qualifier (optional)',
    'GO ID',
    'DB:Reference(|DB:Reference)',
    'Evidence',
    'With (or) From (optional)',
    'Aspect',
    'DB_Object_Name(|Name) (optional)',
    'DB_Object_Synonym(|Synonym) (optional)',
    'DB_Object_Type',
    'taxon(|taxon)',
    'Date',
    'Assigned_by',
    'Annotation_Extension (optional)',
    'Gene_Product_Form_ID (optional)')
)
# mmu_go[,V1:=NULL]
mmu_go[,ID:=`GO ID`]
mmu_go[,UNIPROT:=`DB_Object_ID`]
mmu_go[,SYMBOL:=DB_Object_Symbol]
invisible(mmu_go[Aspect=='P',Aspect:='BP'])
invisible(mmu_go[Aspect=='F',Aspect:='MF'])
invisible(mmu_go[Aspect=='C',Aspect:='CC'])
mmu_go[,ONTOLOGY:=paste0('GO', '_', Aspect)]
mmu_go <- mmu_go[,list(UNIPROT, SYMBOL, ID, ONTOLOGY, present = T)]
invisible(lapply(colnames(mmu_go), function(i) mmu_go[eval(as.name(i))=='',eval(i):=NA]))

setkey(mmu_idmap, UNIPROT)
setkey(mmu_go, UNIPROT)
mmu_go <- merge(mmu_go, mmu_idmap, by='UNIPROT', allow.cartesian = T, all = T)
mmu_go[is.na(SYMBOL.x),SYMBOL.x:=SYMBOL.y]
mmu_go[,SYMBOL.y:=NULL]
setnames(mmu_go, 'SYMBOL.x', 'SYMBOL')
mmu_go <- mmu_go[!is.na(ID),]


mmu_kegg <- fread('dat/annotations/keggedges_mmu', header = F)
setnames(mmu_kegg, c('ID','ENTREZID'))
mmu_kegg[,ID:=gsub('^path\\:', '', ID)]
mmu_kegg[,ENTREZID:=gsub('^mmu\\:', '', ENTREZID)]
setcolorder(mmu_kegg, c('ENTREZID', 'ID'))
invisible(lapply(colnames(mmu_kegg), function(i) mmu_kegg[eval(as.name(i))=='',eval(i):=NA]))

setkey(mmu_idmap, ENTREZID)
setkey(mmu_kegg, ENTREZID)
mmu_kegg <- merge(mmu_kegg, mmu_idmap, by='ENTREZID', allow.cartesian = T, all = T)
# mmu_kegg[is.na(SYMBOL.x),SYMBOL.x:=SYMBOL.y]
# mmu_kegg[,SYMBOL.y:=NULL]
# setnames(mmu_go, 'SYMBOL.x', 'SYMBOL')
mmu_kegg[,ONTOLOGY:='KEGG']
mmu_kegg <- mmu_kegg[!is.na(ID)]
mmu_kegg[,present:=T]


mmu_rct <- fread(
  'dat/annotations/reactomepathways_mmu',
  stringsAsFactors = F,
  integer64 = 'double',
  col.names = c(
    'Source_database_identifier',
    'Reactome_Stable_identifier',
    'URL',
    'Event_(Pathway_or_Reaction)_Name',
    'Evidence_Code',
    'Species'
  )
)

mmu_rct <- mmu_rct[Species=='Mus musculus',]

setnames(mmu_rct, c('Source_database_identifier', 'Reactome_Stable_identifier','Event_(Pathway_or_Reaction)_Name'), c('UNIPROT','ID','Name'))
mmu_rct <- mmu_rct[,list(UNIPROT, ID, Name)]
invisible(lapply(colnames(mmu_rct), function(i) mmu_rct[eval(as.name(i))=='',eval(i):=NA]))

setkey(mmu_idmap, UNIPROT)
setkey(mmu_rct, UNIPROT)
mmu_rct <- merge(mmu_rct, mmu_idmap, by='UNIPROT', allow.cartesian = T, all.x = T)
mmu_rct[,ONTOLOGY:='REACTOME']
mmu_rct[,present:=T]
mmu_rct_names <- mmu_rct[,list(ID, Name)]
setkey(mmu_rct_names)
mmu_rct_names <- unique(mmu_rct_names)
mmu_rct <- mmu_rct[,list(UNIPROT, SYMBOL, ENTREZID, ID, ONTOLOGY, present)]
setkey(mmu_rct)
mmu_rct <- unique(mmu_rct)


mmu <- rbindlist(list(mmu_go, mmu_kegg, mmu_rct), use.names = T)[!is.na(present),]

setkey(mmu)
mmu <- unique(mmu)



# # PREVIOUS VERSION
# # sapply(columns(org.Hs.eg.db), function(i) head(keys(org.Hs.eg.db, keytype = i)))
# mmu <- data.table(select(x = org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = 'UNIPROT'), columns=c('ENTREZID', 'SYMBOL', 'PATH', 'GO'), keytype = 'UNIPROT'))
# mmu[,EVIDENCE:=NULL]
# mmu <- mmu[!is.na(PATH) | !is.na(GO) | !is.na(ONTOLOGY),]
# setkey(mmu)
# mmu <- unique(mmu)
# setnames(mmu, 'GO', 'ID')
# mmu <- rbindlist(list(mmu[!is.na(PATH),list(UNIPROT,ENTREZID,SYMBOL,ID=PATH,ONTOLOGY='KEGG')],mmu[!is.na(ID),list(UNIPROT,ENTREZID,SYMBOL,ID,ONTOLOGY=paste0('GO_', ONTOLOGY))]))
# setkey(mmu)
# mmu <- unique(mmu)
# mmu[ONTOLOGY=='KEGG',ID:=paste0('mmu',ID)]
# mmu[,present:=T]
# # mmu_func <- castfun('UNIPROT',mmu)
# # sapply(columns(reactome.db), function(i) head(keys(reactome.db, keytype = i)))
# rct <- data.table(select(x = reactome.db, keys = keys(reactome.db, keytype = 'ENTREZID'), columns=c('ENTREZID', 'PATHID', 'PATHNAME'), keytype = 'ENTREZID'))
# rct[,org:=unlist(lapply(strsplit(PATHNAME, ': '), function(i) i[[1]]))]
# rct[,PATHNAME:=unlist(lapply(strsplit(PATHNAME, ': '), function(i) i[[2]]))]
# rct <- rct[org%in%c('Homo sapiens'),]
# setnames(rct, colnames(rct), c('ENTREZID', 'ID', 'Name', 'org'))
# rct[,org:=NULL]
# rct[,ONTOLOGY:='Reactome']
# 
# mrg <- merge(rct, mmu[,list(ENTREZID,UNIPROT,SYMBOL)], by='ENTREZID', all.x = T, allow.cartesian = T)
# setkey(mrg)
# mrg <- unique(mrg)
# mrg[,present:=T]
# 
# mmu <- rbindlist(list(mmu, mrg[,list(UNIPROT, ENTREZID, SYMBOL, ID, ONTOLOGY, present)]))


# ---------------------------------------- #
# Names
# ---------------------------------------- #

h <- basicTextGatherer()
curlPerform(url='http://rest.kegg.jp/list/pathway/ath', writefunction = h$update)
kegg_names_ath <- fread(h$value(), header = F)
setnames(kegg_names_ath, colnames(kegg_names_ath), c('ID', 'Name'))
kegg_names_ath[,ID:=gsub('path\\:', '', ID)]
rm(h)
h <- basicTextGatherer()
curlPerform(url='http://rest.kegg.jp/list/pathway/hsa', writefunction = h$update)
kegg_names_hsa <- fread(h$value(), header = F)
setnames(kegg_names_hsa, colnames(kegg_names_hsa), c('ID', 'Name'))
kegg_names_hsa[,ID:=gsub('path\\:', '', ID)]
rm(h)
h <- basicTextGatherer()
curlPerform(url='http://rest.kegg.jp/list/pathway/mmu', writefunction = h$update)
kegg_names_mmu <- fread(h$value(), header = F)
setnames(kegg_names_mmu, colnames(kegg_names_mmu), c('ID', 'Name'))
kegg_names_mmu[,ID:=gsub('path\\:', '', ID)]
rm(h)


GO <- data.table(select(GO.db, keys=keys(GO.db, keytype = 'GOID'), columns=c("GOID","ONTOLOGY", "TERM"), keytype="GOID"))
setnames(GO, colnames(GO), c('ID', 'ONTOLOGY', 'Name'))
GO[,ONTOLOGY:=paste0('GO_', ONTOLOGY)]
fwrite(GO[,list(ID,ONTOLOGY)], 'dat/term_to_ontology.csv', sep = '\t', col.names = T, row.names = F)
GO[,ONTOLOGY:=NULL]
GO <- GO[GOID!='all']

# func_names <- rbindlist(list(kegg_names,GO))

go_offspring <- c(as.list(GOMFOFFSPRING), as.list(GOCCOFFSPRING), as.list(GOBPOFFSPRING))

go_children <- c(as.list(GOMFCHILDREN), as.list(GOCCCHILDREN), as.list(GOBPCHILDREN))


# mrg_names <- mrg[,list(ID, Name)]
# setkey(mrg_names)
# mrg_names <- unique(mrg_names)

# func_names <- rbindlist(list(func_names, ath_rct_names, hsa_rct_names, mmu_rct_names, hsa_cor_names, ath_cor_names))

# GOs <- as.list(GOTERM)

ath_ann <- list(names=rbindlist(list(GO,kegg_names_ath,ath_rct_names,ath_cor_names)), go_offspring=go_offspring, go_children=go_children, annotations=ath, idmap=ath_idmap)
hsa_ann <- list(names=rbindlist(list(GO,kegg_names_hsa,hsa_rct_names,hsa_cor_names)), go_offspring=go_offspring, go_children=go_children, annotations=hsa, idmap=hsa_idmap)
mmu_ann <- list(names=rbindlist(list(GO,kegg_names_mmu,mmu_rct_names)), go_offspring=go_offspring, go_children=go_children, annotations=mmu, idmap=mmu_idmap)

ontomap <- fread('dat/term_to_ontology.csv', stringsAsFactors = F, integer64 = 'double')
ath_ontomap <- unique(ath_ann$annotations[,list(ID,ONTOLOGY)])
hsa_ontomap <- unique(hsa_ann$annotations[,list(ID,ONTOLOGY)])
mmu_ontomap <- unique(mmu_ann$annotations[,list(ID,ONTOLOGY)])

ath_ontomap <- unique(rbindlist(list(ath_ontomap, ontomap)))
hsa_ontomap <- unique(rbindlist(list(hsa_ontomap, ontomap)))
mmu_ontomap <- unique(rbindlist(list(mmu_ontomap, ontomap)))

ath_ann$ontomap <- ath_ontomap
hsa_ann$ontomap <- hsa_ontomap
mmu_ann$ontomap <- mmu_ontomap

saveRDS(ath_ann, 'preproc/ath_ann.rds', compress = T)
saveRDS(hsa_ann, 'preproc/hsa_ann.rds', compress = T)
saveRDS(mmu_ann, 'preproc/mmu_ann.rds', compress = T)

# saveRDS(golist, 'preproc/golist.rds', compress = T)

