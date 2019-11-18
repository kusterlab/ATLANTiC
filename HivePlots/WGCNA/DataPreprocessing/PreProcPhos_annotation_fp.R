require(data.table)
annotation <- fread('/media/general/projects/NCI_60_phospho/Final_dataset/MQ1.5.5.1_CRC65_FP_complete/txt/proteinGroups.txt', stringsAsFactors = F, integer64 = 'double')
invisible(lapply(colnames(annotation), function(i) annotation[eval(as.name(i))=='',eval(i):='NA']))

annotation[,firstprot:=unlist(
  lapply(
    strsplit(
      `Majority protein IDs`,
      ';'),
    function(i)
      head(i,
           1)
  )
)
]

annotation[,lengthprot:=unlist(
  lapply(
    strsplit(
      `Majority protein IDs`,
      ';'),
    function(i)
      length(i)
  )
)
]

annotation[,firstgene:=unlist(
  lapply(
    strsplit(
      `Gene names`,
      ';'),
    function(i)
      head(i,
           1)
  )
)
]

annotation[,lengthgene:=unlist(
  lapply(
    strsplit(
      `Gene names`,
      ';'),
    function(i)
      length(i)
  )
)
]


# annotation[,firstpos:=unlist(
#   lapply(
#     strsplit(
#       Positions,
#       ';'),
#     function(i)
#       head(i,
#            1)
#   )
# )
# ]

# annotation[,lengthpos:=unlist(
#   lapply(
#     strsplit(
#       Positions,
#       ';'),
#     function(i)
#       length(i)
#   )
# )
# ]

# annotation[,firstposgene:=unlist(
#   lapply(
#     strsplit(
#       `Positions within proteins`,
#       ';'),
#     function(i)
#       head(i,
#            1)
#   )
# )
# ]

# annotation[,lengthposgene:=unlist(
#   lapply(
#     strsplit(
#       `Positions within proteins`,
#       ';'),
#     function(i)
#       length(i)
#   )
# )
# ]

# annotation[firstgene=='NA',firstgene:='NA']
annotation[`Gene names`=='NA',`Gene names`:='']
# annotation[firstpos=='NA',firstpos:='']

annotation[,rid:=paste0(
  id,
  ' ',
  `Gene names`
)]


# annotation[,label:=paste0(firstgene, '_p', `Amino acid`, firstpos)]
setkey(annotation, rid)
saveRDS(annotation, '/media/msdata5/users_files/martin/postdoc/projects/phosphoproject/dat/preproc/crc65_fp_annotation.rds', compress = T)
