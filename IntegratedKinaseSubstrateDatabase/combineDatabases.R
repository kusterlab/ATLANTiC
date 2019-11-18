pnet <- read.delim("Res/20160803_kinaseSubstrateDatabaseParse/phosphoNetwork.selftab.txt", stringsAsFactors = FALSE)
pelm <- read.delim("Res/20160803_kinaseSubstrateDatabaseParse/phosphoELM.selftab.txt", stringsAsFactors = FALSE)
psp <- read.delim("Res/20160803_kinaseSubstrateDatabaseParse/phosphoSitePlus.selftab.txt", stringsAsFactors = FALSE)
signor <- read.delim("Res/20160803_kinaseSubstrateDatabaseParse/signor.selftab.txt", stringsAsFactors = FALSE)
unip <- read.delim("Res/20160803_kinaseSubstrateDatabaseParse/uniprot.selftab.txt", stringsAsFactors = FALSE)
kea <- read.delim("Res/20160803_kinaseSubstrateDatabaseParse/KEA.selftab.txt", stringsAsFactors = FALSE)


head(kea)

integ <- data.frame(kinase = c(pnet$kinase,
                      pelm$kinases,
                      psp$GENE,
                      signor$EntityA,
                      unip$kinase,
                      kea$X
                      ),
           sub_acc = c(sapply(strsplit(pnet$ACC, ";"), "[", 1),
                       pelm$acc,
                       sapply(strsplit(psp$SUB_ACC_ID, "-"), "[", 1),
                       signor$IdB,
                       unip$proteinAcc,
                       kea$acc1
                       ),
           sub_geneName = c(pnet$protein, 
                            pelm$geneName,
                            psp$SUB_GENE,
                            signor$EntityB,
                            unip$geneName,
                            sapply(strsplit(kea$x, "_"), "[", 1)
                            ),
           sub_uniprotID = c(pnet$ID, 
                             pelm$ID,
                             psp$id,
                             signor$ID,
                             unip$protein,
                             kea$ID
                             ),
           sub_seq7 = c(pnet$flankseq7, 
                        substr(pelm$sequence, 9, 23),
                        psp$SITE_...7_AA,
                        signor$MechanismSequences,
                        substr(unip$flankseq15, 9, 23),
                        substr(kea$flankseq15, 9, 23)
                        ),
           aminoAcid = c(pnet$aminoAcid, 
                         pelm$code,
                         psp$aminoAcid,
                         signor$aminoAcid,
                         unip$aminoAcid,
                         kea$aa
                         ),
           position = c(pnet$position, 
                        pelm$position,
                        substr(psp$SUB_MOD_RSD, 2, nchar(psp$SUB_MOD_RSD)),
                        signor$position,
                        unip$site,
                        kea$pos
                        ),
           database = c(rep("phosphoNetwork", nrow(pnet)), 
                        rep("phosphoELM", nrow(pelm)),
                        rep("phosphoSitePlus", nrow(psp)),
                        rep("signor", nrow(signor)),
                        rep("uniprot", nrow(unip)),
                        rep("KEA", nrow(kea))
                        )
           )


integ$aminoAcid <- toupper(integ$aminoAcid)
integ$sub_acc <- as.character(integ$sub_acc)
integ$sub_uniprotID <- as.character(integ$sub_uniprotID)
integ$sub_geneName <- toupper(integ$sub_geneName)
integ$sub_seq7 <- toupper(integ$sub_seq7)
integ$position <- as.integer(as.character(integ$position))
integ$kinase <- toupper(integ$kinase)
integ$kinase[which(integ$kinase == "AUTOCATALYSIS")] <- integ$sub_geneName[which(integ$kinase == "AUTOCATALYSIS")]

sapply(integ, class)
integ <- integ[substr(integ$sub_seq7, 8, 8) == integ$aminoAcid &
                 integ$aminoAcid %in% c("S", "T", "Y"),  ]

dim(integ)
write.table(integ, file="Res/20160803_kinaseSubstrateDatabaseParse/integDatabase.txt",
            col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")
saveRDS(integ, "Res/20160803_kinaseSubstrateDatabaseParse/integDatabase.RDS")
length(unique(integ$kinase))


tail(sort(table(paste(integ$kinase, integ$sub_seq7))), n=100)
length(unique(paste(integ$kinase, integ$sub_seq7)))














