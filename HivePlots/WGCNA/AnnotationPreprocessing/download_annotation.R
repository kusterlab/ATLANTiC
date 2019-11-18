require(RCurl)
require(R.utils)
require(data.table)

# URLs
urls <- fread('dat/todownload.csv', stringsAsFactors = F, integer64 = 'double')

links <- urls$URL
zipfiles <- sapply(links, function(i) tail(strsplit(i, split = '\\/')[[1]], 1))
urls[,Name2:=paste0(Name, ifelse(lapply(strsplit(zipfiles, '\\.'), tail, 1)=='gz','.gz', ''))]
urls[,Name2:=paste0(Name, ifelse(lapply(strsplit(zipfiles, '\\.'), tail, 1)=='zip','.zip', ''))]
# files <- gsub('.gz', '', zipfiles)
# names(files) <- urls$Name
names(zipfiles) <- urls$Name
names(links) <- urls$Name
files <- urls$Name2
names(files) <- urls$Name
downloaded <- list.files(file.path('dat', 'annotations'), recursive = F, pattern = paste0(paste0(paste0('^', names(files)), '$'), collapse = '|'))
notdownloaded <- names(files)[!names(files)%in%downloaded]

# Download & Unzip
if(length(notdownloaded)>0) {
  for (i in notdownloaded) {
    f = CFILE(file.path('dat', 'annotations', files[i]), mode="wb")
    curlPerform(url = links[i], writedata = f@ref)
    close(f)
    try(gunzip(filename = file.path('dat', 'annotations', files[i])), silent = T)
    try({
        unzip(zipfile = file.path('dat', 'annotations', files[i]), exdir = file.path('dat', 'annotations'))
        file.remove(file.path('dat', 'annotations', files[i]))
      }, silent = T)
  }
}
