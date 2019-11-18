al2hex <- function(alphacol, bg = 'white') {
  compcol <- col2rgb(alphacol, alpha = T)
  compbg <- col2rgb(bg, alpha = F)
  compcol <- rbind(compcol[4,]*compcol[1:3,,drop=F]/255, 255-compcol[4,,drop = F])
  rescol <- compcol[1:3,,drop=F] + (compcol[4,]*compbg)/255
  return(rgb(t(rescol), maxColorValue = 255))
}

library(extrafont)
# font_import()
# fonts()
# fonttable()
loadfonts(device="postscript")
loadfonts(device="pdf")

.cs <- list()
.cs <- within(.cs, {
  colorScheme <- function(scheme = c("c2_1", "c2_2", "c3", "c4", "c5", "c6"), names = NULL) {
    tumcolors <- c(TUMmauve = "#69085A",
                   TUMviolet = "#0F1B5F",
                   TUMblackblue = "#003359",
                   TUMdarkblue = "#005293",
                   TUMblue = "#0073CF",
                   TUMlightblue = "#64A0C8",
                   TUMciel = "#98C6EA",
                   TUMturk = "#00778A",
                   TUMgreen = "#007C30",
                   TUMlightgreen = "#679A1D",
                   TUMolive = "#A2AD00",
                   TUMbeige = "#DAD7CB",
                   TUMyellow = "#FFDC00",
                   TUMgold = "#F9BA00",
                   TUMorange = "#E37222",
                   TUMbrick = "#D64C13",
                   TUMred = "#C4071B",
                   TUMblood = "#9C0D16")
    c2_1 <- c("TUMblue", "TUMred")
    c2_2 <- c("TUMturk", "TUMbeige")
    c3 <- c("TUMmauve", "TUMgreen", "TUMbrick")
    c4 <- c("TUMmauve", "TUMlightblue", "TUMlightgreen", "TUMbrick")
    c5 <- c("TUMmauve", "TUMdarkblue", "TUMciel", "TUMlightgreen", "TUMgold")
    c6 <- c(c5, "TUMred")
    
    col <- tumcolors[get(scheme)]
    if (!is.null(names))
      names(col) <- names
    col
  }
  
  
  
  
  
  
  crcColor <- rgb(227, 114, 34, maxColorValue = 255 )
  crcColor <- structure(
    rep(crcColor, 65), 
    names = c("C10_CRC65", "C106_CRC65", "C125-PM_CRC65", "C32_CRC65", "C70_CRC65", "C75_CRC65", "C80_CRC65", 
              "C84_CRC65", "C99_CRC65", "CaCo-2_CRC65", "CaR-1_CRC65", "CC07_CRC65", "CC20_CRC65", 
              "CCK-81_CRC65", "CoCM-1_CRC65", "Colo 320DM_CRC65", "Colo 741_CRC65", "Colo-678_CRC65", 
              "DLD-1_CRC65", "GP2d_CRC65", "HCA-46_CRC65", "HCA-7_CRC65", "HCT116_CRC65", "HDC-111_CRC65", 
              "HDC-114_CRC65", "HDC-135_CRC65", "HDC-142_CRC65", "HDC-143_CRC65", "HDC-54_CRC65", "HDC-57_CRC65", 
              "HDC-73_CRC65", "HDC-8_CRC65", "HDC-82_CRC65", "HDC-9_CRC65", "HRA-19_CRC65", "HT29_CRC65", 
              "HT55_CRC65", "LIM1863_CRC65", "LoVo_CRC65", "LS 174T_CRC65", "LS 180_CRC65", "LS1034_CRC65", 
              "LS123_CRC65", "LS411_CRC65", "LS513_CRC65", "NCI-H548_CRC65", "NCI-H716_CRC65", "NCI-H747_CRC65", 
              "OXCO-1_CRC65", "OXCO-3_CRC65", "PC/JW_CRC65", "RCM-1_CRC65", "RKO_CRC65", "SK-CO-1_CRC65", 
              "SNU-C2B_CRC65", "SW1116_CRC65", "SW1417_CRC65", "SW403_CRC65", "SW48_CRC65", "SW480_CRC65", 
              "SW620_CRC65", "SW837_CRC65", "SW948_CRC65", "T84_CRC65", "VACO 4A_CRC65"))
  
  tooColor <- c("RE" = "magenta4","LC" = "gray50","CO" = "orange3","PR" = "black",
                "OV" = "red3","LE" = "cyan3", "CN" = "green4","ME" = "brown","BR" = "blue3")
  nciColor <- c("786O_NCI60" = "magenta4",
                A498_NCI60 = "magenta4",
                A549_NCI60 = "gray50",
                ACHN_NCI60 = "magenta4",
                BT549_NCI60 ="blue3",
                CAKI1_NCI60 ="magenta4",
                CCRFCEM_NCI60 = "cyan3",
                COLO205_NCI60 = "orange3",
                DU145_NCI60 ="black",
                EKVX_NCI60 = "gray50",
                HCC2998_NCI60 = "orange3",
                HCT116_NCI60 = "orange3",
                HCT15_NCI60 = "orange3",
                HL60_NCI60 = "cyan3",
                HOP62_NCI60 = "gray50",
                HOP92_NCI60 = "gray50",
                HS578T_NCI60 = "blue3",
                HT29_NCI60 = "orange3",
                IGROV1_NCI60 = "red3",
                K562_NCI60 = "cyan3",
                KM12_NCI60 = "orange3",
                LOXIMVI_NCI60 = "brown",
                M14_NCI60 = "brown",
                MALME3M_NCI60 = "brown",
                MCF7_NCI60 = "blue3",
                MDAMB231_NCI60 = "blue3",
                MDAMB435_NCI60 = "brown",
                MDAMB468_NCI60 = "blue3",
                MOLT4_NCI60 = "cyan3",
                NCIADRES_NCI60 = "red3",
                NCIH226_NCI60 = "gray50",
                NCIH23_NCI60 = "gray50",
                NCIH322M_NCI60 = "gray50",
                NCIH460_NCI60 = "gray50",
                NCIH522_NCI60 = "gray50",
                OVCAR3_NCI60 = "red3",
                OVCAR4_NCI60 = "red3",
                OVCAR5_NCI60 = "red3",
                OVCAR8_NCI60 = "red3",
                PC3_NCI60 = "black",
                RPMI8226_NCI60 = "cyan3",
                RXF393_NCI60 = "magenta4",
                SF268_NCI60 = "green4",
                SF295_NCI60 = "green4",
                SF539_NCI60 = "green4",
                SKMEL2_NCI60 = "brown",
                SKMEL28_NCI60 = "brown",
                SKMEL5_NCI60 = "brown",
                SKOV3_NCI60 = "red3",
                SN12C_NCI60 = "magenta4",
                SNB19_NCI60 = "green4",
                SNB75_NCI60 = "green4",
                SR_NCI60 = "cyan3",
                SW620_NCI60 = "orange3",
                T47D_NCI60 = "blue3",
                TK10_NCI60 = "magenta4",
                U031_NCI60 = "magenta4",
                U251_NCI60 = "green4",
                UACC257_NCI60 = "brown",
                UACC62_NCI60 = "brown")
  
  tumblue <- "#0073CF"
  tColor <- function(x, alpha = 100) {
    rgbs <- col2rgb(x)
    vv <- rgb(rgbs[1, ], rgbs[2, ], rgbs[3, ], alpha = alpha, maxColorValue = 255)
    names(vv) <- names(x)
    vv
  }
})
