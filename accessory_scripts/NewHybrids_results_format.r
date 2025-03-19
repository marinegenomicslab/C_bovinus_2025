#!/bin/R

## Usage ##
#NewHybrids_results_format.r <Group 1 name> <File of Group 1 indvs> <Group 2 name> <File of Group 2 indvs>

#Importing data
args <- commandArgs(trailingOnly=TRUE)
dat <- read.table("aa-ScaledLikelihood_indv.txt", head=T, sep="\t")
GRP1 <- as.character(as.matrix(read.table(args[2], head=F)[,1]))
GRP2 <- as.character(as.matrix(read.table(args[4], head=F)[,1]))
GRP1.sum <- apply(dat[dat$Indv %in% GRP1, 3:7], 2, sum)
GRP2.sum <- apply(dat[dat$Indv %in% GRP2, 3:7], 2, sum)
names.list <- colnames(dat)

colnames(dat)[which(GRP1.sum == max(GRP1.sum))+2] <- args[1]
colnames(dat)[which(GRP2.sum == max(GRP2.sum))+2] <- args[3]
colnames(dat)[colnames(dat) == "X0.000.0.500.0.500.0.000"] <- "F1_cross"
colnames(dat)[colnames(dat) == "X0.000.0.250.0.250.0.500"] <- paste(colnames(dat)[names.list == "X0.000.0.000.0.000.1.000"],"BX")
colnames(dat)[colnames(dat) == "X0.500.0.250.0.250.0.000"] <- paste(colnames(dat)[names.list == "X1.000.0.000.0.000.0.000"],"BX")

write.table(dat[, -2], file="aa-ScaledLikelihood_edit.txt", quote=F, col.names=T, row.names=F, sep="\t")
### End of NewHybrids_results_format.r ###
