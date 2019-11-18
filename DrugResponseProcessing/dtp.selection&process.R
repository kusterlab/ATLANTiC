dtp <- list()
dtp <- within (dtp, {
  rdir1 <- "/media/kusterlab/users_files/Chen Meng/Projects/PhosphoNCI60/"
  rdir2 <- "/media/kusterlab/internal_projects/active/NCI60_Phospho/Chen/"
  
  dtp <- read.delim(file.path(rdir1, "Data/DTP_DrugSensitivity/DTP_NCI60.xls/DTP_NCI60_1_5_1.txt"), 
                    comment.char = "#", check.names = FALSE, as.is = TRUE)
  kbs <- read.delim(file.path(rdir2, "Dat/kbscreen&review/kinobeadScreen.txt"))
  santosRev <- read.delim(file.path(rdir2, "Dat/kbscreen&review/santos_nrd.2016.230-s2.txt"))
  
  
  ###
  drug <- list(santosRev = santosRev$PARENT_PREF_NAME, 
               kinoScreen = kbs$Drug,
               dtp = dtp$Name)
  
  drug <- lapply(drug, as.character)
  drug <- lapply(drug, unique)
  drug.l <- lapply(drug, tolower)
  drug.nop <- lapply(drug.l, function(x) gsub("[[:punct:]]| ", "", x))
  
  map.santosRev <- dtp$Name %in% drug$dtp[drug.nop$dtp %in% drug.nop$santosRev]
  map.kinoScreen <- dtp$Name %in% drug$dtp[drug.nop$dtp %in% drug.nop$kinoScreen]
  
  ic <- c(1:4, 65:66)
  sel <- map.santosRev | map.kinoScreen | !dtp$`FDA status` %in% c("", "-")
  dtp.sel <- cbind(dtp[sel, ic], santosRev = as.integer(map.santosRev[sel]), 
                   kinoScreen = as.integer(map.kinoScreen[sel]), dtp[sel, -ic])
  
  dtp.annot <- dtp.sel[, 1:8]
  dtp.gi50 <- apply(dtp.sel[, -(1:8)], 2, as.numeric)
  
  rownames(dtp.annot) <- rownames(dtp.gi50) <- paste(dtp$NSC[sel], dtp$Name[sel], sep = "_")
  rm(list = c("dtp.sel", "sel", "ic", "map.kinoScreen", "map.santosRev", "drug.nop",
              "drug.l", "drug", "santosRev", "kbs", "dtp"))
})
