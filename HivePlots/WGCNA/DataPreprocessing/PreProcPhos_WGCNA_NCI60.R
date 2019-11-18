source('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/src/WGCNA.R')
options(expressions = 1e5)

allsites_nci60_tryp <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_nci60_tryp/RDS/allsites_incnas.rds')
allsites_nci60_tryp_red <- allsites_nci60_tryp[rowSums(!is.na(allsites_nci60_tryp))>=7,]

func_ann <- readRDS('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/preproc/hsa_ann.rds')
traits_nci60 <- readRDS('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/phosphoproject/dat/traits_nci60.rds')

doWGCNA(mat1 = allsites_nci60_tryp_red,
        mat2 = NULL,
        setLabels = c('NoImputation'),
        type = 'signed',
        traits = traits_nci60,
        func_ann = func_ann,
        identifier = 'SYMBOL',
        softPower = 12, #3
        MEDissThres = 0.1,
        agglomeration = 'average',
        colord = 'sorted',
        selectedTraits = colnames(traits_nci60),
        plottedTraits = colnames(traits_nci60)[1:9],
        outpath = '/media/kusterlab/users_files/Martin Frejno/postdoc/projects/phosphoproject/res/20180315_wgcna_allsites_allproteins_incnas_nci60',
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
        nodesize = 10,
        cutoff = 0.05,
        gsubid = '_p[STY][0-9]*(_[0-9])?$')

rm(allsites_nci60_tryp)
rm(allsites_nci60_tryp_red)
gc()