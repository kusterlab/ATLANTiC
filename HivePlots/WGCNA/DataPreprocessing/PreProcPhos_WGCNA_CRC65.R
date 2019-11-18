source('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/src/WGCNA.R')
options(expressions = 1e5)

allsites_crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180109_wgcna_allsites_allproteins_crc65/RDS/allsites_incnas.rds')
allsites_crc65_red <- allsites_crc65[rowSums(!is.na(allsites_crc65))>=7,]

func_ann <- readRDS('/media/kusterlab/users_files/Martin Frejno/postdoc/projects/arabidopsis/preproc/hsa_ann.rds')
traits_crc65 <- readRDS('/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/dat/traits_crc65.rds')

doWGCNA(mat1 = allsites_crc65_red,
        mat2 = NULL,
        setLabels = c('NoImputation'),
        type = 'signed',
        traits = traits_crc65,
        func_ann = func_ann,
        identifier = 'SYMBOL',
        softPower = 12, #3
        MEDissThres = 0.1,
        agglomeration = 'average',
        colord = 'sorted',
        selectedTraits = colnames(traits),
        plottedTraits = colnames(traits)[1:26],
        outpath = '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180315_wgcna_allsites_allproteins_incnas_crc65/',
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

rm(allsites_crc65)
rm(allsites_crc65_red)
gc()