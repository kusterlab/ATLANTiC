cli <- readRDS("Res/20170324_Fluororacil_TCGA/cli.drug.RDS")
cap <- unique(grep("cape|Xeloda", cli$all6cancer$drug_name, ignore.case = TRUE, value = TRUE))
capecitabine <- cli$all6cancer[cli$all6cancer$drug_name %in% cap,] 
fu5 <- cli$fu5

fu5$drug <- "5FU"
capecitabine$drug <- "capecitabine"
doi <- rbind(fu5, capecitabine)


table(doi$measure_of_response)
difres <- list()
difres$CR <- doi[doi$measure_of_response == "Complete Response", ]
difres$CPD <- doi[doi$measure_of_response == "Clinical Progressive Disease", ]
difres$PR <- doi[doi$measure_of_response == "Partial Response", ]
difres$SD <- doi[doi$measure_of_response == "Stable Disease", ]


dff <- do.call("rbind", difres)
ct <- split(dff$bcr_patient_barcode, dff$CancerType)
saveRDS(dff, file = "Res/20170324_Fluororacil_TCGA/tcga.patient.drug.RDS")

######

library(MASS)
y <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
x1 <- y + rnorm(9, sd = 0.2)
x2 <- x1[c(4:6, 1:3, 7:9)]
plot(x1, y)
plot(x2, y)

m <- polr(as.factor(y) ~ x1)
summary(m)
ctable <- coef(summary(m))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2


library(RTCGAToolbox)


cbind(ctable, "p value" = p)







