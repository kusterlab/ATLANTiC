require(data.table)
require(reshape2)

# source('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/src/MQPSMMap.R')
# genPSMtable('/media/kusterlab/internal_projects/active/NCI60_Phospho/Final_dataset/MQ1.5.5.1_pNCI60_pCRC65_complete/txt/')
# genPSMtable('/media/kusterlab/internal_projects/active/NCI60_Phospho/Final_dataset/MQ1.5.5.1_CRC65_FP_complete/txt/')
# genPSMtable('/media/kusterlab/internal_projects/active/NCI60_Phospho/Final_dataset/MQ1.5.5.1_NCI60_FP_complete/txt/')

phosprots <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/Final_dataset/MQ1.5.5.1_pNCI60_pCRC65_complete/txt/proteinGroups.txt', integer64 = 'double')
setkey(phosprots, Reverse)
phosprots <- phosprots[!.("+")]
setkey(phosprots, `Only identified by site`)
phosprots <- phosprots[!.("+")]
invisible(phosprots[,id:=as.character(id)])
phosprotsgn <- lapply(strsplit(phosprots$`Gene names`, ';'), head, 1)
phosprotsgnlen <- lapply(phosprotsgn, length)
for (i in seq_along(phosprotsgnlen)) {
  if (phosprotsgnlen[[i]]==0) {
    phosprotsgn[[i]] <- NA
  }
}
phosprots[,firstgene:=unlist(phosprotsgn)]
phosprots[`Gene names`=='NA',`Gene names`:='']
phosprots[,rid:=paste0(id,' ',`Gene names`)]

fwrite(phosprots, '~/phospho/Manuscript/current/supplementTables/new/phosprots_annotation.txt', sep = '\t', col.names = T, row.names = F)

prots_crc65 <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/Final_dataset/MQ1.5.5.1_CRC65_FP_complete/txt/proteinGroups.txt', integer64 = 'double')
setkey(prots_crc65, Reverse)
prots_crc65 <- prots_crc65[!.("+")]
setkey(prots_crc65, `Only identified by site`)
prots_crc65 <- prots_crc65[!.("+")]
invisible(prots_crc65[,id:=as.character(id)])
prots_crc65gn <- lapply(strsplit(prots_crc65$`Gene names`, ';'), head, 1)
prots_crc65gnlen <- lapply(prots_crc65gn, length)
for (i in seq_along(prots_crc65gnlen)) {
  if (prots_crc65gnlen[[i]]==0) {
    prots_crc65gn[[i]] <- NA
  }
}
prots_crc65[,firstgene:=unlist(prots_crc65gn)]
raw_crc65 <- as.matrix(prots_crc65[,grep('Intensity ', colnames(prots_crc65)),with=F])
rownames(raw_crc65) <- prots_crc65$id
rawdt_crc65 <- data.table(melt(raw_crc65, varnames = c('PROTID', 'EXPID'), value.name = 'Raw Intensity'))
rawdt_crc65[,PROTID:=as.character(PROTID)]
rawdt_crc65[,EXPID:=gsub('Intensity ', '', as.character(EXPID))]
setkey(rawdt_crc65, EXPID, PROTID)


prots_nci60 <- fread('/media/kusterlab/internal_projects/active/NCI60_Phospho/Final_dataset/MQ1.5.5.1_NCI60_FP_complete/txt/proteinGroups.txt', integer64 = 'double')
setkey(prots_nci60, Reverse)
prots_nci60 <- prots_nci60[!.("+")]
setkey(prots_nci60, `Only identified by site`)
prots_nci60 <- prots_nci60[!.("+")]
invisible(prots_nci60[,id:=as.character(id)])
prots_nci60gn <- lapply(strsplit(prots_nci60$`Gene names`, ';'), head, 1)
prots_nci60gnlen <- lapply(prots_nci60gn, length)
for (i in seq_along(prots_nci60gnlen)) {
  if (prots_nci60gnlen[[i]]==0) {
    prots_nci60gn[[i]] <- NA
  }
}
prots_nci60[,firstgene:=unlist(prots_nci60gn)]
raw_nci60 <- as.matrix(prots_nci60[,grep('Intensity ', colnames(prots_nci60)),with=F])
rownames(raw_nci60) <- prots_nci60$id
rawdt_nci60 <- data.table(melt(raw_nci60, varnames = c('PROTID', 'EXPID'), value.name = 'Raw Intensity'))
rawdt_nci60[,PROTID:=as.character(PROTID)]
rawdt_nci60[,EXPID:=gsub('Intensity ', '', as.character(EXPID))]
setkey(rawdt_nci60, EXPID, PROTID)



mqpsmmap <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Final_dataset/MQ1.5.5.1_pNCI60_pCRC65_complete/txt/evmap_red.rds')
setkey(mqpsmmap, EXPID, PROTID)

ev_phosprotein <- mqpsmmap[`Phospho (STY)`>0,list(Intensity=sum(Intensity, na.rm = T)),by=list(EXPID,PROTID)]
setkey(ev_phosprotein, PROTID)
setkey(phosprots, id)
ev_phosprotein[,`Gene names`:=phosprots[PROTID,firstgene]]

ev_phosprotein_crc65 <- ev_phosprotein[grepl('crc65', EXPID, ignore.case = T),]

ev_phosprots_export_crc65 <- acast(data = ev_phosprotein_crc65, formula = PROTID~EXPID, value.var = 'Intensity')
ev_phosprots_export_crc65 <- ev_phosprots_export_crc65[as.character(sort(as.numeric(rownames(ev_phosprots_export_crc65)))),]
setkey(phosprots, id)
rownames(ev_phosprots_export_crc65) <- phosprots[.(rownames(ev_phosprots_export_crc65)),rid]
colnames(ev_phosprots_export_crc65) <- gsub('_CRC65', '', colnames(ev_phosprots_export_crc65))
ev_phosprots_export_crc65 <- ev_phosprots_export_crc65[,sort(colnames(ev_phosprots_export_crc65))]
ev_phosprots_export_crc65[ev_phosprots_export_crc65==0] <- NA
ev_phosprots_export_crc65 <- log2(ev_phosprots_export_crc65)
ev_phosprots_export_crc65 <- data.table(ID=rownames(ev_phosprots_export_crc65), ev_phosprots_export_crc65)
fwrite(ev_phosprots_export_crc65, '~/phospho/Manuscript/current/supplementTables/new/phosprots_abundance_crc65.txt', sep = '\t', col.names = T, row.names = F)

phosprotid_crc65 <- ev_phosprotein_crc65[,mean(Intensity),by=list(PROTID,`Gene names`)][,.SD[which.max(V1),PROTID],by=`Gene names`][,V1]
ev_phosprotein_crc65 <- ev_phosprotein_crc65[PROTID%in%phosprotid_crc65,]
ev_phosprots_crc65 <- acast(data = ev_phosprotein_crc65, formula = `Gene names`~EXPID, value.var = 'Intensity')
ev_phosprots_crc65[ev_phosprots_crc65==0] <- NA
ev_phosprots_crc65 <- log2(ev_phosprots_crc65)


ev_phosprotein_nci60 <- ev_phosprotein[grepl('nci60', EXPID, ignore.case = T),]

ev_phosprots_export_nci60 <- acast(data = ev_phosprotein_nci60, formula = PROTID~EXPID, value.var = 'Intensity')
ev_phosprots_export_nci60 <- ev_phosprots_export_nci60[as.character(sort(as.numeric(rownames(ev_phosprots_export_nci60)))),]
setkey(phosprots, id)
rownames(ev_phosprots_export_nci60) <- phosprots[.(rownames(ev_phosprots_export_nci60)),rid]
colnames(ev_phosprots_export_nci60) <- gsub('_NCI60', '', colnames(ev_phosprots_export_nci60))
ev_phosprots_export_nci60 <- ev_phosprots_export_nci60[,sort(colnames(ev_phosprots_export_nci60))]
ev_phosprots_export_nci60[ev_phosprots_export_nci60==0] <- NA
ev_phosprots_export_nci60 <- log2(ev_phosprots_export_nci60)
ev_phosprots_export_nci60 <- data.table(ID=rownames(ev_phosprots_export_nci60), ev_phosprots_export_nci60)
fwrite(ev_phosprots_export_nci60, '~/phospho/Manuscript/current/supplementTables/new/phosprots_abundance_nci60.txt', sep = '\t', col.names = T, row.names = F)

phosprotid_nci60 <- ev_phosprotein_nci60[,mean(Intensity),by=list(PROTID,`Gene names`)][,.SD[which.max(V1),PROTID],by=`Gene names`][,V1]
ev_phosprotein_nci60 <- ev_phosprotein_nci60[PROTID%in%phosprotid_nci60,]
ev_phosprots_nci60 <- acast(data = ev_phosprotein_nci60, formula = `Gene names`~EXPID, value.var = 'Intensity')
ev_phosprots_nci60[ev_phosprots_nci60==0] <- NA
ev_phosprots_nci60 <- log2(ev_phosprots_nci60)

mqpsmmap_crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Final_dataset/MQ1.5.5.1_CRC65_FP_complete/txt/evmap_red.rds')
setkey(mqpsmmap_crc65, EXPID, PROTID)

ev_protein_crc65 <- mqpsmmap_crc65[,list(Intensity=sum(Intensity, na.rm = T)),by=list(EXPID,PROTID)]
setkey(ev_protein_crc65, PROTID)
setkey(prots_crc65, id)
ev_protein_crc65[,`Gene names`:=prots_crc65[PROTID,firstgene]]
protid_crc65 <- ev_protein_crc65[,mean(Intensity),by=list(PROTID,`Gene names`)][,.SD[which.max(V1),PROTID],by=`Gene names`][,V1]
ev_protein_crc65 <- ev_protein_crc65[PROTID%in%protid_crc65,]

ev_prots_crc65 <- acast(data = ev_protein_crc65[grepl('CRC65', EXPID),], formula = `Gene names`~EXPID, value.var = 'Intensity')
ev_prots_crc65[ev_prots_crc65==0] <- NA
ev_prots_crc65 <- log2(ev_prots_crc65)

setkey(ev_protein_crc65, EXPID, PROTID)
setkey(rawdt_crc65, EXPID, PROTID)
merged_crc65 <- ev_protein_crc65[rawdt_crc65,]
merged_crc65[`Raw Intensity`==0,`Raw Intensity`:=NA]
merged_crc65[,Intensity:=signif(Intensity, digits = 5)]
merged_crc65[,cor(`Raw Intensity`,Intensity, use = 'p')]


mqpsmmap_nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/Final_dataset/MQ1.5.5.1_NCI60_FP_complete/txt/evmap_red.rds')
setkey(mqpsmmap_nci60, EXPID, PROTID)

ev_protein_nci60 <- mqpsmmap_nci60[,list(Intensity=sum(Intensity, na.rm = T)),by=list(EXPID,PROTID)]
setkey(ev_protein_nci60, PROTID)
prots_nci60[,id:=as.character(id)]
setkey(prots_nci60, id)
ev_protein_nci60[,`Gene names`:=prots_nci60[PROTID,firstgene]]
protid_nci60 <- ev_protein_nci60[,mean(Intensity),by=list(PROTID,`Gene names`)][,.SD[which.max(V1),PROTID],by=`Gene names`][,V1]
ev_protein_nci60 <- ev_protein_nci60[PROTID%in%protid_nci60,]

ev_prots_nci60 <- acast(data = ev_protein_nci60[grepl('NCI60', EXPID),], formula = `Gene names`~EXPID, value.var = 'Intensity')
ev_prots_nci60[ev_prots_nci60==0] <- NA
ev_prots_nci60 <- log2(ev_prots_nci60)

setkey(ev_protein_nci60, EXPID, PROTID)
setkey(rawdt_nci60, EXPID, PROTID)
merged_nci60 <- ev_protein_nci60[rawdt_nci60,]
merged_nci60[`Raw Intensity`==0,`Raw Intensity`:=NA]
merged_nci60[,Intensity:=signif(Intensity, digits = 5)]
merged_nci60[,cor(`Raw Intensity`,Intensity, use = 'p')]

# protid_crc65 <- intersect(rownames(ev_prots_crc65), rownames(ev_phosprots_crc65))
cellid_crc65 <- intersect(colnames(ev_prots_crc65), colnames(ev_phosprots_crc65))
# protid_nci60 <- intersect(rownames(ev_prots_nci60), rownames(ev_phosprots_nci60))
cellid_nci60 <- intersect(colnames(ev_prots_nci60), colnames(ev_phosprots_nci60))

ev_prots_crc65 <- ev_prots_crc65[,cellid_crc65]
ev_phosprots_crc65 <- ev_phosprots_crc65[,cellid_crc65]
ev_prots_nci60 <- ev_prots_nci60[,cellid_nci60]
ev_phosprots_nci60 <- ev_phosprots_nci60[,cellid_nci60]

saveRDS(ev_prots_crc65, file = '/media/kusterlab/users_files/Martin Frejno/postdoc/projects/phosphoproject/dat/preproc/ev_prots_crc65.rds', compress = T)
saveRDS(ev_phosprots_crc65, file = '/media/kusterlab/users_files/Martin Frejno/postdoc/projects/phosphoproject/dat/preproc/ev_phosprots_crc65.rds', compress = T)
saveRDS(ev_prots_nci60, file = '/media/kusterlab/users_files/Martin Frejno/postdoc/projects/phosphoproject/dat/preproc/ev_prots_nci60.rds', compress = T)
saveRDS(ev_phosprots_nci60, file = '/media/kusterlab/users_files/Martin Frejno/postdoc/projects/phosphoproject/dat/preproc/ev_phosprots_nci60.rds', compress = T)
