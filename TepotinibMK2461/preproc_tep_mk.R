require(tidyverse)
require(tidyselect)

give_names <- function(.) {
  return(c('HDC-8_1', 'HDC-8_2', 'HDC-8_3'))
}

give_conc <- function(.) {
  return('Dose [uM]')
}

basedir <- '/media/kusterlab/internal_projects/active/NCI60_Phospho/martin/res/20180413_tepotinib_mk/'
files <- list.files(basedir, pattern = '.txt', full.names = F)

dosedata <- do.call(bind_rows, lapply(files, function(i) {
  path <- file.path(basedir, i)
  path %>%
    read_tsv() %>%
    rename_if(.predicate = grepl('\\[uM\\]', colnames(.)), .funs = give_conc) %>%
    rename_if(.predicate = grepl('HDC-8', colnames(.)), .funs = give_names) %>%
    gather(key = 'Cell line', value = 'Viability', contains('HDC-8')) %>%
    mutate(Replicate = unlist(lapply(str_split(string = `Cell line`, pattern = '_'), tail, 1))) %>%
    mutate(`Cell line` = unlist(lapply(str_split(string = `Cell line`, pattern = '_'), head, 1))) %>%
    mutate(Treatment = switch(gsub('.txt', '', i),
                              combination_mkscale = 'Combination (MK-2461-scale)',
                              combination_tepscale = 'Combination (Tepotinib-scale)',
                              mk2461_mkscale = 'MK-2461',
                              tepotinib_tepscale = 'Tepotinib')) %>%
    dplyr::select(Treatment, `Dose [uM]`, `Cell line`, Replicate, Viability)
}))

write_tsv(dosedata, path = file.path(basedir, 'dosedata.tsv'))
