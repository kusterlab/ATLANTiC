slist <- readRDS("Res/20170110_mq15processLOD2/sites.trypsin.RDS")

library(stringr)
library(pheatmap)

###
dis <- read.delim("Dat/phosphoSitePlus/20180315/Disease-associated_sites", comment.char = "#", stringsAsFactors = FALSE, header = TRUE)
sites <- read.delim("Dat/phosphoSitePlus/20180315/Phosphorylation_site_dataset", comment.char = "#", stringsAsFactors = FALSE)


### antibody
pppid <- paste(str_split_fixed(sites$ACC_ID, "-|\\.", 2)[, 1], toupper(sites$SITE_...7_AA), sep = "_")
withab <- pppid[sites$CST_CAT != "" & sites$ORGANISM == "human"]

### disease associated
disid <- paste(str_split_fixed(dis$ACC_ID, "-|\\.", 2)[, 1], toupper(dis$SITE_...7_AA), sep = "_")

### measured ID
mid <- paste(str_split_fixed(slist$annot$`Leading proteins`, "-|\\.|;", 2)[, 1], 
             substr(slist$annot$`Sequence window`, start = 9, stop = 23), sep = "_")

################################# get matrices ####################################
i_disease <- rownames(slist$annot)[mid %in% disid]
i_antibody <- rownames(slist$annot)[mid %in% withab]
# i_regulate <- rownames(slist$annot)[slist$annot$`Regulatory site` == "+"]


mat_nci60 <- slist$fun_prepElnet(x = slist, panel = "NCI60")$imputed
mat_crc65 <- slist$fun_prepElnet(x = slist, panel = "crc65", celllines = setdiff(colnames(slist$crc65_expr), "CoCM-1_CRC65"))$imputed

mat_dis_nci60 <- mat_nci60[rownames(mat_nci60) %in% i_disease, ]
mat_dis_crc65 <- mat_crc65[rownames(mat_crc65) %in% i_disease, ]

mat_ab_nci60 <- mat_nci60[rownames(mat_nci60) %in% i_antibody, ]
mat_ab_crc65 <- mat_crc65[rownames(mat_crc65) %in% i_antibody, ]

mat_reg_nci60 <- slist$fun_prepElnet(x = slist, panel = "NCI60", regsitesOnly = TRUE)$imputed
mat_reg_crc65 <- slist$fun_prepElnet(x = slist, panel = "crc65", regsitesOnly = TRUE, 
                                     celllines = setdiff(colnames(slist$crc65_expr), "CoCM-1_CRC65"))$imputed


######################## heatmap 1 ############################ 
s <- unique(slist$nci60_xref[, c("too", "color")])
parlist <- list(scale = "row", 
                clustering_distance_rows = "correlation", 
                clustering_distance_cols = "correlation", 
                clustering_method = "ward.D",
                # breaks = c(-10, seq(-3, 3, length.out = 99), 10), 
                cellwidth = 10, cellheight = 0.5,
                annotation_col = data.frame(too = slist$nci60_xref[colnames(mat_ab_nci60), "too"], 
                                            row.names = colnames(mat_ab_nci60)), 
                annotation_colors = list(too = structure(s$color, names = s$too)))


#
ip <- parlist
ip$mat <- mat_reg_nci60
ip$filename = "Res/20180316_wgcnaFigures/heatmap_reg_nci60.tiff"
do.call(pheatmap, ip)

ip <- parlist
ip$mat <- mat_ab_nci60
ip$filename = "Res/20180316_wgcnaFigures/heatmap_ab_nci60.tiff"
do.call(pheatmap, ip)

ip <- parlist
ip$mat <- mat_dis_nci60
ip$filename = "Res/20180316_wgcnaFigures/heatmap_dis_nci60.tiff"
do.call(pheatmap, ip)

#
ip <- parlist
ip$mat <- mat_reg_crc65
ip$filename = "Res/20180316_wgcnaFigures/heatmap_reg_crc65.tiff"
do.call(pheatmap, ip)

ip <- parlist
ip$mat <- mat_ab_crc65
ip$filename = "Res/20180316_wgcnaFigures/heatmap_ab_crc65.tiff"
do.call(pheatmap, ip)

ip <- parlist
ip$mat <- mat_dis_crc65
ip$filename = "Res/20180316_wgcnaFigures/heatmap_dis_crc65.tiff"
do.call(pheatmap, ip)




######################## heatmap 2 ############################ 
s <- unique(slist$nci60_xref[, c("too", "color")])
parlist <- list(scale = "row", 
                clustering_distance_rows = "correlation", 
                clustering_distance_cols = "correlation", 
                clustering_method = "ward.D",
                breaks = c(-10, seq(-3, 3, length.out = 99), 10), 
                cellwidth = 10, cellheight = 0.5,
                annotation_col = data.frame(too = slist$nci60_xref[colnames(mat_ab_nci60), "too"], 
                                            row.names = colnames(mat_ab_nci60)), 
                annotation_colors = list(too = structure(s$color, names = s$too)))


#
ip <- parlist
ip$mat <- mat_reg_nci60
ip$filename = "Res/20180316_wgcnaFigures/heatmap_reg_nci60_break.tiff"
do.call(pheatmap, ip)

ip <- parlist
ip$mat <- mat_ab_nci60
ip$filename = "Res/20180316_wgcnaFigures/heatmap_ab_nci60_break.tiff"
do.call(pheatmap, ip)

ip <- parlist
ip$mat <- mat_dis_nci60
ip$filename = "Res/20180316_wgcnaFigures/heatmap_dis_nci60_break.tiff"
do.call(pheatmap, ip)

#
ip <- parlist
ip$mat <- mat_reg_crc65
ip$filename = "Res/20180316_wgcnaFigures/heatmap_reg_crc65_break.tiff"
do.call(pheatmap, ip)

ip <- parlist
ip$mat <- mat_ab_crc65
ip$filename = "Res/20180316_wgcnaFigures/heatmap_ab_crc65_break.tiff"
do.call(pheatmap, ip)

ip <- parlist
ip$mat <- mat_dis_crc65
ip$filename = "Res/20180316_wgcnaFigures/heatmap_dis_crc65_break.tiff"
do.call(pheatmap, ip)

