source('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/src/WGCNA.R')
options(expressions = 1e5)

ev_prots_nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/ev_prots_nci60.rds')
ev_phosprots_nci60 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/preproc/ev_phosprots_nci60.rds')
ids <- names(which(rowSums(!is.na(ev_prots_nci60))>=7 & rowSums(!is.na(ev_phosprots_nci60))>=7))
ev_prots_nci60_red <- ev_prots_nci60[ids,]
ev_phosprots_nci60_red <-ev_phosprots_nci60[ids,]
identical(dimnames(ev_prots_nci60_red),dimnames(ev_phosprots_nci60_red))

func_ann <- readRDS('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/preproc/hsa_ann.rds')
traits_nci60 <- readRDS('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/phosphoproject/dat/traits_nci60.rds')

doWGCNA(mat1 = ev_prots_nci60_red,
        mat2 = ev_phosprots_nci60_red,
        setLabels = c('Proteins', 'Phosphoproteins'),
        type = 'signed',
        traits = traits_nci60,
        func_ann = func_ann,
        identifier = 'SYMBOL',
        softPower = 10, #3
        MEDissThres = 0.1,
        agglomeration = 'average',
        deepSplit = 4,
        pamRespectsDendro = T,
        colord = 'sorted',
        selectedTraits = colnames(traits),
        plottedTraits = colnames(traits)[1:26],
        outpath = '/media/kusterlab/users_files/Martin Frejno/postdoc/projects/phosphoproject/res/20180405_wgcna_prot_phosprot_incnas_nci60/',
        scale = T,
        scaleP = 0.5,
        scalecut = 0,
        sampleord = 'hclust',
        excludeoutliers = F,
        removenas = F,
        minNOBS = 7,
        filtermad = F,
        quantcut = 0,
        heatmaps = F,
        writeentirenetwork = F,
        automatedtomcut = F,
        tomcut = 0.02,
        minModuleSize = 30,
        cutoff = 0.05)

rm(ev_prots_nci60)
rm(ev_prots_nci60_red)
rm(ev_phosprots_nci60)
rm(ev_phosprots_nci60_red)
gc()