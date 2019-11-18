library(openxlsx)
mm <- read.xlsx("/media/kusterlab/users_files/Runsheng Zheng/1_RTK project/2_Results/12_Competition_pull-down/1_Peptide_list_for_synthesis/tab_filtered_2_with Biogrid and PSP (Autosaved).xlsx")

tab <- read.delim("Res/20180328_outlierAnalysis/outlierTable.txt", stringsAsFactors = FALSE)

#####################
grep("ALK", unique(mm$pep), value = TRUE)
pd <- mm$Gene.names[mm$pep == "iBAQ_2315_B01_P024581_R1_ALK_1507"]

ot <- tab$Gene.name[tab$panel == "CRC65" & tab$cell == "C10"]
ot <- sapply(strsplit(ot, ";"), "[", 1)

intersect(pd, ot)



#####################
grep("ALK", unique(mm$pep), value = TRUE)
pd <- mm$Gene.names[mm$pep == "iBAQ_2315_B01_P024581_R1_ALK_1507"]

ot <- tab$Gene.name[tab$panel == "CRC65" & tab$cell == "C10"]
ot <- sapply(strsplit(ot, ";"), "[", 1)

intersect(pd, ot)


#####################
grep("EGFR", unique(mm$pep), value = TRUE)
pd <- mm$Gene.names[mm$pep == "iBAQ_1886_A03_P019636_R1_EGFR_998"]

ot <- tab$Gene.name[tab$panel == "CRC65" & tab$cell == "HDC-8"]
ot <- sapply(strsplit(ot, ";"), "[", 1)

intersect(pd, ot)


#####################
grep("DDR1", unique(mm$pep), value = TRUE)
pd <- mm$Gene.names[mm$pep %in% c("iBAQ_2297_H02_P024561_R1_DDR1_797", 
                                  "iBAQ_1632_A12_P017068_R1_DDR1_520")]

ot <- tab$Gene.name[tab$panel == "CRC65" & tab$cell == "HDC-8"]
ot <- sapply(strsplit(ot, ";"), "[", 1)

intersect(pd, ot)


#####################
grep("FGFR2", unique(mm$pep), value = TRUE)
pd <- mm$Gene.names[mm$pep %in% c("iBAQ_1768_F05_P018527_R1_FGFR2_566")]

ot <- tab$Gene.name[tab$panel == "CRC65" & tab$cell == "NCI-H716"]
ot <- sapply(strsplit(ot, ";"), "[", 1)

intersect(pd, ot)
