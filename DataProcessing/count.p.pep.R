library(data.table)
###
tryp <- fread("../Final_dataset/MQ1.5.5.1_pNCI60_pCRC65_complete/txt/evidence.txt", integer64 = "double", data.table = FALSE)
gluc <- fread("../Final_dataset/MQ1.5.5.1_pNCI60_GluC_complete/txt/evidence.txt", integer64 = "double", data.table = FALSE)

###
tryp.seq <- gsub("\\(ac\\)|\\(ox\\)", "", grep("\\(ph\\)", tryp$`Modified sequence`, value = TRUE))
gluc.seq <- gsub("\\(ac\\)|\\(ox\\)", "", grep("\\(ph\\)", gluc$`Modified sequence`, value = TRUE))

length(unique(c(tryp.seq, gluc.seq)))


##
a <- tryp$`Modified sequence`[grepl("\\(ph\\)", tryp$`Modified sequence`)  & tryp$Score > 40]
a <- gsub("\\(ac\\)|\\(ox\\)", "", a)
b <- gluc$`Modified sequence`[grepl("\\(ph\\)", gluc$`Modified sequence`)  & gluc$Score > 40]
b <- gsub("\\(ac\\)|\\(ox\\)", "", b)

length(unique(a))
length(unique(b))
length(unique(c(a, b)))


ll <- intersect(a, b)
length(ll)
head(a)
head(b)












