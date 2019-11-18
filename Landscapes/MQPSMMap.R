# pathtotxtfolder <- '/media/kusterlab/internal_projects/active/Arabidopsis_proteome/Arabdiopsis thaliana/Full proteome/Data_Tissues/combined_MBR_30T_FPandPhospho/'
# pathtomqpar <- '/media/kusterlab/internal_projects/active/Arabidopsis_proteome/Arabdiopsis thaliana/Maxquant Search AT_FP_P/mqpar.xml'

genPSMtable <- function(pathtotxtfolder) { # , pathtomqpar
  require(data.table)
  require(parallel)
  # require(xml2)
  # require(dplyr)
  require(reshape2)
  alldup <- function (value) { 
    duplicated(value) | duplicated(value, fromLast = TRUE)
  }
  
  # varmods <- unique(xml_find_all(read_xml(pathtomqpar), "//parameterGroups/parameterGroup/variableModifications/string") %>% xml_text())
  # quantmods <- unique(xml_find_all(read_xml(pathtomqpar), "//restrictMods/string") %>% xml_text())
  # excmods <- setdiff(varmods, quantmods)
  
  cat('Read in proteinGroups!\n')
  proteinGroups <- fread(file.path(pathtotxtfolder, 'proteinGroups.txt'), integer64 = 'double')
  setkey(proteinGroups, Reverse)
  proteinGroups <- proteinGroups[!.("+")]
  setkey(proteinGroups, `Only identified by site`)
  proteinGroups <- proteinGroups[!.("+")]
  invisible(proteinGroups[,id:=as.character(id)])
  
  # proteinGroups[,firstid:=unlist(lapply(strsplit(`Protein IDs`, ';'), function(i) i[1]))]
  
  
  # Arabidopsis-specific
  # proteinGroups[,`Gene names`:=unlist(lapply(strsplit(`Majority protein IDs`, ';'), function(i) paste(unique(gsub('\\.[0-9]*$', '', i)), collapse = ';')))]
  setkey(proteinGroups, id)
  
  cat('Read in peptides!\n')
  peptides <- fread(file.path(pathtotxtfolder, 'peptides.txt'), integer64 = 'double')
  invisible(peptides[,id:=as.character(id)])
  setkey(peptides, Reverse)
  peptides <- peptides[!.("+")]
  setkey(peptides, `Unique (Groups)`)
  peptides <- peptides[.('yes'),]
  
  # for (i in excmods) {
  #   peptides <- peptides[!grepl(pattern = "[0-9]+", x = eval(as.symbol(paste0(i, ' site IDs')))),]
  # }
  
  cat('Find razor peptide assignment!\n')
  pepraz <- proteinGroups[,list(id, `Peptide IDs`,`Peptide is razor`)]
  pepids <- strsplit(pepraz[,`Peptide IDs`], ';')
  idx <- strsplit(pepraz[,`Peptide is razor`], ';')
  idx <- lapply(idx, function(i) as.logical(toupper(i)))
  pepids_red <- mclapply(seq_along(pepids), function(i) pepids[[i]][idx[[i]]], mc.cores = 10)
  names(pepids_red) <- pepraz$id
  pepmap <- data.table(PROTID=unlist(mclapply(names(pepids_red), function(i) rep(i, times=length(pepids_red[[i]])), mc.cores = 10)),PEPID=unlist(pepids_red))
  setkey(pepmap, PEPID)
  # pepmap <- pepmap[.(peptides$id),]
  
  cat('Read in evidence!\n')
  evidence <- fread(file.path(pathtotxtfolder, 'evidence.txt'), integer64 = 'double')
  invisible(evidence[,id:=as.character(id)])
  invisible(evidence[,`Peptide ID`:=as.character(`Peptide ID`)])
  setkey(evidence, Reverse)
  evidence <- evidence[!J("+")]
  mods <- evidence[,unique(gsub('[0-9] ', '', unlist(strsplit(Modifications, ','))))]
  mods <- mods[!mods%in%c("Unmodified")]
  setkey(evidence, `Peptide ID`)
  evidence <- evidence[.(pepmap[,PEPID]),]
  evidence <- evidence[!is.na(Modifications),]
  setkey(evidence, `Peptide ID`)
  
  # ppeps <- evidence[`Phospho (STY)`>0,unique(`Peptide ID`)]
  # evidence <- evidence[!.(ppeps),]
  
  # for (i in excmods) {
  #   evidence <- evidence[eval(as.symbol(i))==0,]
  # }
  
  # evidence <- evidence[,list(id, `Peptide ID`, `Retention time`, `Retention time calibration`, `Calibrated retention time`, `m/z`, Charge, `Uncalibrated mass error [ppm]`, Score, Experiment, Fraction, Intensity)]
  # # evidence <- evidence[,list(id,Reverse,`Potential contaminant`,`Modified sequence`,Proteins, `Leading proteins`, `Leading razor protein`,`Protein group IDs`, `Peptide ID`, `Gene names`,Experiment,Fraction,Intensity)]
  # gc()
  
  colstoget <- c(RT='Retention time', CRT='Calibrated retention time', RTC='Retention time calibration', MZ='m/z', CHARGE='Charge', MerrPPM='Uncalibrated mass error [ppm]', SCORE='Score', EVID='id', PEPID='Peptide ID', EXPID='Experiment', Intensity='Intensity')
  colspres <- unlist(sapply(colstoget, function(i) grep(gsub("\\]", "\\\\]", gsub("\\[", "\\\\[", gsub("\\/", "\\\\/", paste0('^', i, '$')))), colnames(evidence), ignore.case = T, value = T)))
  if (!identical(colstoget,colspres)) {
    warning('evidence.txt has unexpected column names!')
  }
  
  cat('Generate evidence PSM Map!\n')
  setnames(evidence, colspres, names(colspres))
  evmap <- evidence[,list(RT, CRT, RTC, MZ, CHARGE, MerrPPM, SCORE, EVID, PEPID, EXPID, Intensity),nomatch=0]
  evmap[,c(mods):=data.table(evidence[,sapply(mods, function(i) eval(as.name(i)))])]
  evmap[,Unmodified:=as.integer(rowSums(sapply(mods, function(i) eval(as.name(i))))==0)]
  # evmap <- evidence[pepmap[,PEPID],list(EVID=id, PEPID=`Peptide ID`, EXPID=Experiment, REV=Reverse, CON=`Potential contaminant`, Intensity=Intensity),nomatch=0]
  setkey(pepmap, PEPID)
  # pepmap <- pepmap[!.(ppeps),]
  setkey(evmap, PEPID)
  evmap <- pepmap[evmap,]
  # invisible(evmap[,PROTID:=pepmap[PEPID,PROTID]])
  
  # setkey(proteinGroups, id)
  # invisible(evmap[,GENEID:=proteinGroups[PROTID,`Gene names`]])
  
  # setkey(evmap, EXPID)
  # 
  digits <- 5
  # invisible(sapply(grep("Intensity ", colnames(proteinGroups), value=T), function(i) proteinGroups[eval(as.name(i))==0,c(i):=NA]))
  # samp <- evmap[,unique(EXPID)[30]]
  # setkey(proteinGroups, id)
  # plot(log2(evmap[.(samp),][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(Intensity, na.rm=T), digits = digits),by=PROTID][,V1]),log2(proteinGroups[evmap[.(samp),][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(Intensity, na.rm=T), digits = digits),by=PROTID][,PROTID],eval(as.symbol(paste0('Intensity ', samp)))]), pch='.')

  # setkey(evmap, GENEID)
  # evmap <- evmap[!J(""),]
  setkey(evmap, EXPID)
  
  evmap_red <- evmap[!is.na(Intensity),]
  # evmap_red_cptac <- evmap_red[unique(EXPID)[nchar(unique(EXPID))==2],]
  # evmap_red_crc64 <- evmap_red[unique(EXPID)[nchar(unique(EXPID))!=2],]
  
  # gc()

  # require(Biostrings)
  # fasta <- readAAStringSet(filepath = pathtofasta, format = 'fasta', use.names = T)
  # fst <- data.table(prt=unlist(lapply(strsplit(names(fasta), '\\|'), function(i) i[2])),
  #                   gn=unlist(lapply(strsplit(unlist(lapply(strsplit(names(fasta), 'GN='), function(i) i[2])), ' '), function(j) j[1])))
  
  cat('Generate raw intensity matrix for filtering!\n')
  rawc <- proteinGroups[,c("id", grep("Intensity ", colnames(proteinGroups), value = T)),with=F]
  setnames(rawc, colnames(rawc), gsub("Intensity ", "", colnames(rawc)))
  rawcm <- as.matrix(rawc[,!colnames(rawc)%in%c("id"),with=F])
  rownames(rawcm) <- rawc[,`id`]
  rawcm[rawcm==0] <- NA
  rawcm <- rawcm[which(rowSums(is.na(rawcm))!=ncol(rawcm)),]
  setkey(evmap_red, PROTID)
  evmap_red[.(rownames(rawcm)),NOPGINT:=F]
  evmap_red[is.na(NOPGINT),NOPGINT:=T]
  
  if (any(grepl("iBAQ ", colnames(proteinGroups)))) {
    cat('Generate iBAQ intensity matrix for filtering!\n')
    ibaqc <- proteinGroups[,c("id", grep("iBAQ ", colnames(proteinGroups), value = T)),with=F]
    setnames(ibaqc, colnames(ibaqc), gsub("iBAQ ", "", colnames(ibaqc)))
    ibaqcm <- as.matrix(ibaqc[,!colnames(ibaqc)%in%c("id"),with=F])
    rownames(ibaqcm) <- ibaqc[,`id`]
    ibaqcm[ibaqcm==0] <- NA
    rawc <- proteinGroups[,c("id", grep("Intensity ", colnames(proteinGroups), value = T)),with=F]
    setnames(rawc, colnames(rawc), gsub("Intensity ", "", colnames(rawc)))
    rawcm <- as.matrix(rawc[,!colnames(rawc)%in%c("id"),with=F])
    rownames(rawcm) <- rawc[,`id`]
    rawcm[rawcm==0] <- NA
    cat('Calculate expepted number of peptides!\n')
    expeplist <- apply(round(rawcm/ibaqcm), 1, function(i) unique(i)[!is.na(unique(i))])
    rawcm <- rawcm[!rownames(rawcm)%in%names(which(unlist(lapply(expeplist, length)==0))),]
    ibaqcm <- ibaqcm[!rownames(ibaqcm)%in%names(which(unlist(lapply(expeplist, length)==0))),]
    expeplist_red <- expeplist[!names(expeplist)%in%names(which(unlist(lapply(expeplist, length)==0)))]
    proteinGroups_red <- proteinGroups[!.(names(which(unlist(lapply(expeplist, length)==0))))]
    expep <- data.table(#`Gene names`=proteinGroups[,`Gene names`],
      id=proteinGroups_red[,id],
      # accs=proteinGroups_red[,firstid],
      peps=unlist(expeplist_red))
    
    # setkey(evmap_red, GENEID)
    # setkey(expep_red, `Gene names`)
    # invisible(evmap_red[,iBAQ_gn:=Intensity/expep_red[GENEID,peps]])
    setkey(evmap_red, PROTID)
    # evmap_red <- evmap_red[.(rownames(rawcm)),]
    evmap_red[.(rownames(rawcm)),EXPEP:=T]
    evmap_red[is.na(EXPEP),EXPEP:=F]
    setkey(expep, id)
    invisible(evmap_red[,iBAQ_pr:=Intensity/expep[PROTID,peps]])
    
    invisible(evmap_red[,compid:=paste(EXPID, PROTID, sep='_|_')])
    setkey(evmap_red, compid)
    evmap_prot <- evmap_red[,signif(sum(Intensity, na.rm=T), digits = digits), by=compid]
    setkey(evmap_prot, compid)
    invisible(evmap_red[,Intfrac_prot:=Intensity/evmap_prot[compid,V1]])
    
    invisible(evmap_red[,compid:=NULL])
    
    setkey(evmap_red, EXPID)
    
    # plot(log2(evmap_red[samp,][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(Intensity, na.rm=T), digits = digits),by=PROTID][,V1]),log2(rawcm[evmap_red[samp,][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(Intensity, na.rm=T), digits = digits),by=PROTID][,PROTID],samp]), pch='.')
    # plot(log2(evmap_red[samp,][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(iBAQ_pr, na.rm=T), digits = digits),by=PROTID][,V1]),log2(ibaqcm[evmap_red[samp,][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(Intensity, na.rm=T), digits = digits),by=PROTID][,PROTID],samp]), pch='.')
  }
  
  cat('Assign razor peptide ids!\n')
  # IMPORTANT
  setkey(evmap_red, PEPID)
  evmap_red[.(unique(peptides$id)),ISRAZ:='F']
  evmap_red[is.na(ISRAZ),ISRAZ:='T']
  setkey(evmap_red, EXPID)
  # plot(log2(evmap_red[ISRAZ=='F',][.(samp),][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(Intensity, na.rm=T), digits = digits),by=PROTID][,V1]),log2(rawcm[evmap_red[ISRAZ=='F',][.(samp),][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(Intensity, na.rm=T), digits = digits),by=PROTID][,PROTID],samp]), pch='.')
  # plot(log2(evmap_red[ISRAZ=='F',][.(samp),][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(iBAQ_pr, na.rm=T), digits = digits),by=PROTID][,V1]),log2(ibaqcm[evmap_red[ISRAZ=='F',][.(samp),][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(Intensity, na.rm=T), digits = digits),by=PROTID][,PROTID],samp]), pch='.')
  
  # rawmatrix <- acast(formula = PROTID~EXPID, data = evmap_red[!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(Intensity, na.rm=T), digits = digits),by=list(EXPID,PROTID)], value.var = 'V1')
  # rawmatrix <- log2(rawmatrix[sort(rownames(rawmatrix)),sort(colnames(rawmatrix))])
  # 
  # genes <- intersect(rownames(rawmatrix),rownames(rawcm))
  # # sanity check
  # diag(cor(rawmatrix[genes,],log2(rawcm[genes,]),use = 'p'))
  # samplemedians <- apply(rawmatrix, 2, median, na.rm=T)
  # overallmedian <- median(rawmatrix, na.rm = T)
  # rawmatrix <- sweep(rawmatrix, 2, samplemedians, `-`)+overallmedian
  # 
  # norazormatrix <- acast(formula = PROTID~EXPID, data = evmap_red[ISRAZ=='F',][!PEPID%in%evmap[`Phospho (STY)`>0,unique(PEPID)],signif(sum(Intensity, na.rm=T), digits = digits),by=list(EXPID,PROTID)], value.var = 'V1')
  # norazormatrix <- log2(norazormatrix[sort(rownames(norazormatrix)),sort(colnames(norazormatrix))])
  # norazormatrix <- sweep(norazormatrix, 2, samplemedians, `-`)+overallmedian
  # saveRDS(norazormatrix, file = 'dat/fp_norazor.rds', compress = T)
  cat('Save PSMmap to disk!\n')
  saveRDS(evmap_red, file = file.path(pathtotxtfolder, 'evmap_red.rds'), compress = T)
}
