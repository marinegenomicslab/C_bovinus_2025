#!/bin/R

## Usage ##
#NewHybrids_vcf_to_txt.r <Input_file> <Export_file>
#input file: vcf file to subset and convert
#export file: txt file of the data foro NewHybrids

#Load libraries
suppressMessages(library(adegenet, quietly=T))
suppressMessages(library(vcfR, quietly=T))

#Load function
Loci_names <- function(NAMES, SEP="[.]", REMOVE=1){
COL <- length(strsplit(head(NAMES,n=1), SEP)[[1]])
TMP_U <- unique(matrix(unlist(strsplit(NAMES,SEP)),ncol=COL,byrow=T)[,1:(COL-REMOVE)])
if(is.matrix(TMP_U)){TMP_DF <- data.frame(TMP_U)
} else {TMP_DF <- data.frame(matrix(TMP_U, ncol=(COL-REMOVE), byrow=T))}
return(tidyr::unite(TMP_DF, "loci", 1:ncol(TMP_DF), sep=SEP))
}

#Load data
args <- commandArgs(trailingOnly=TRUE)
vcf <- read.vcfR(file = args[1], verbose=F)
vcf.gen <- vcfR2genind(vcf)
rm(vcf)
individuals <- indNames(vcf.gen)
locnames <- locNames(vcf.gen)

#Selecting SNPs
keeploc <- sample(locnames, size = 300, replace = FALSE, prob = NULL)
#Making sure there is only one SNP per contig
loc.keep <- Loci_names(keeploc, SEP="_", REMOVE=1)
keeploc <- keeploc[head(as.numeric(rownames(unique(loc.keep))), n=150)]
#Subset vcf file
gen_subset <- vcf.gen[loc = keeploc]

# Use gen_subset to run NewHyrbids
# extract genotype matrix
df_snp <- genind2df(gen_subset, usepop = FALSE)

#missing data is coded for by using 0, so need to update other values
df_snp[df_snp=="11"] <- "22"
df_snp[df_snp=="00"] <- "11"
df_snp[df_snp=="01"] <- "12"
df_snp[df_snp=="10"] <- "21"

# change NA to 0
df_snp[is.na(df_snp)] = 0

# Name individuals numerically
df_snp$LocusNames <- c(1:length(individuals))
#Note: these aren't actually locus names, they are the individuals, but in the data file "Locus Names" $

#Check number of colums in df
print(paste("Number of loci is", ncol(df_snp)-1))

# move sample numbers to first column
df_snp <- df_snp[,c(151, 1:150)]
#write your NewHybrids txt input file
write.table(df_snp, args[2], quote = FALSE, row.names = FALSE)
### End of NewHybrids_vcf_to_txt.r ###