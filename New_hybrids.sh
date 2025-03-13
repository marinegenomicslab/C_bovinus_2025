### Hybrid analysis  of the Gambusia ###
## Earth ##
#Making executable cmd "NewHybrids_vcf_to_txt.r"
{```{R}```
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
#args <- c("~/Workspace/Robert/Gnobilis/analysis/hybrid/Ggei_New_Mexico_indv.recode.vcf", "newhybrids_input_GxN_NM.txt")
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
#150

# move sample numbers to first column
df_snp <- df_snp[,c(151, 1:150)]
#write your NewHybrids txt input file
write.table(df_snp, args[2], quote = FALSE, row.names = FALSE)
}} #notepad cleanup
#Making executable cmd "NewHybrids_results_format.r"
{```{R}```
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
}


## NewHybrids ##
#Taking files with individuals listed by group and combining them for pairwise analysis
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid
#Get data
cp -s ../Sep2023/*.txt .
cp -s ../Sep2023/*.vcf .

#Make affinis and nobilis pairs
cat Gaff.txt Diamond_Y_nobilis.txt > Ga_DY.txt
cat Gaff.txt West_Texas_nobilis.txt > Ga_WTX.txt
cat Gaff.txt New_Mexico_nobilis.txt > Ga_NM.txt

#Make geiseri and nobilis pairs
cat Ggei.txt Diamond_Y_nobilis.txt > Gg_DY.txt
cat Ggei.txt West_Texas_nobilis.txt > Gg_WTX.txt
cat Ggei.txt New_Mexico_nobilis.txt > Gg_NM.txt

#Make affinis and nobilis pairs
cat Gaff.txt Ggei.txt > Ga_Gg.txt

#Make vcf files
vcftools --vcf SNP.TRS.F06.vcf --out Ga_DY --recode --recode-INFO-all --keep Ga_DY.txt
vcftools --vcf SNP.TRS.F06.vcf --out Ga_WTX --recode --recode-INFO-all --keep Ga_WTX.txt
vcftools --vcf SNP.TRS.F06.vcf --out Ga_NM --recode --recode-INFO-all --keep Ga_NM.txt

vcftools --vcf SNP.TRS.F06.vcf --out Gg_DY --recode --recode-INFO-all --keep Gg_DY.txt
vcftools --vcf SNP.TRS.F06.vcf --out Gg_WTX --recode --recode-INFO-all --keep Gg_WTX.txt
vcftools --vcf SNP.TRS.F06.vcf --out Gg_NM --recode --recode-INFO-all --keep Gg_NM.txt

vcftools --vcf SNP.TRS.F06.vcf --out Ga_Gg --recode --recode-INFO-all --keep Ga_Gg.txt

#Adding the characterizations for the analysis
nano ~/bin/genotype_frequencies.txt
{ #Pasted info
   5
 1_Bx      0.00000    0.250000    0.250000   0.50000
 0_Bx      0.50000    0.250000    0.250000   0.00000
 F1        0.00000     .5          .5        0.00000
 Pure_1    0.00000    0.00000     0.00000    1.00000
 Pure_0    1.00000    0.00000     0.00000    0.00000
}
}
## G. geseri vs New Mexico G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid
#Preparing dir
mkdir -p geiseri/NM
cd geiseri/NM
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/Robert/Gambusia/analysis/hybrid/Gg_NM.recode.vcf newhybrids_input_GxN_NM.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Gg_NM.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_GxN_NM.txt;
done

INDV=$(vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Gg_NM.recode.vcf | wc -l);sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" newhybrids_input_GxN_NM.txt;

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_GxN_NM.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui &
sleep 5
cd ..
done
}
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Gg_NM.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.geiseri ~/Workspace/Robert/Gambusia/analysis/hybrid/Ggei.txt G.nobilis ~/Workspace/Robert/Gambusia/analysis/hybrid/New_Mexico_nobilis.txt
cd ..
done
}
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "GxN_NM.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "GxN_NM.xlsx", overwrite = TRUE)
}
}
## G. geseri vs Diamond Y G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid
#Preparing dir
mkdir -p geiseri/DY
cd geiseri/DY
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/Robert/Gambusia/analysis/hybrid/Gg_DY.recode.vcf newhybrids_input_GxN_DY.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Gg_DY.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_GxN_DY.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_GxN_DY.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui &
sleep 5
cd ..
done
}
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Gg_DY.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.geiseri ~/Workspace/Robert/Gambusia/analysis/hybrid/Ggei.txt G.nobilis ~/Workspace/Robert/Gambusia/analysis/hybrid/Diamond_Y_nobilis.txt
cd ..
done
}
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "GxN_DY.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "GxN_DY.xlsx", overwrite = TRUE)
}
}
## G. geseri vs West Texas G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid
#Preparing dir
mkdir -p geiseri/WTX
cd geiseri/WTX
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/Robert/Gambusia/analysis/hybrid/Gg_WTX.recode.vcf newhybrids_input_GxN_WTX.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Gg_WTX.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_GxN_WTX.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_GxN_WTX.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui &
sleep 5
cd ..
done
}
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Gg_WTX.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.geiseri ~/Workspace/Robert/Gambusia/analysis/hybrid/Ggei.txt G.nobilis ~/Workspace/Robert/Gambusia/analysis/hybrid/West_Texas_nobilis.txt
cd ..
done
}
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "GxN_WTX.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "GxN_WTX.xlsx", overwrite = TRUE)
}
}

## G. affinis vs New Mexico G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid
#Preparing dir
mkdir -p affinis/NM
cd affinis/NM
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/Robert/Gambusia/analysis/hybrid/Ga_NM.recode.vcf newhybrids_input_AxN_NM.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Ga_NM.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_AxN_NM.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_AxN_NM.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui &
sleep 5
cd ..
done
}
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Ga_NM.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.affinis ~/Workspace/Robert/Gambusia/analysis/hybrid/Gaff.txt G.nobilis ~/Workspace/Robert/Gambusia/analysis/hybrid/New_Mexico_nobilis.txt
cd ..
done
}
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "AxN_NM.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "AxN_NM.xlsx", overwrite = TRUE)
}
}
## G. affinis vs Diamond Y G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid
#Preparing dir
mkdir -p affinis/DY
cd affinis/DY
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/Robert/Gambusia/analysis/hybrid/Ga_DY.recode.vcf newhybrids_input_AxN_DY.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Ga_DY.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_AxN_DY.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_AxN_DY.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui &
sleep 5
cd ..
done
}
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Ga_DY.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.affinis ~/Workspace/Robert/Gambusia/analysis/hybrid/Gaff.txt G.nobilis ~/Workspace/Robert/Gambusia/analysis/hybrid/Diamond_Y_nobilis.txt
cd ..
done
}
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "AxN_DY.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "AxN_DY.xlsx", overwrite = TRUE)
}
}
## G. geseri vs West Texas G. nobilis ##
{```{bash}```
#Generating data
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid
#Preparing dir
mkdir -p affinis/WTX
cd affinis/WTX
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/Robert/Gambusia/analysis/hybrid/Ga_WTX.recode.vcf newhybrids_input_AxN_WTX.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Ga_WTX.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_AxN_WTX.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_AxN_WTX.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui &
sleep 5
cd ..
done
}
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Ga_WTX.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.affinis ~/Workspace/Robert/Gambusia/analysis/hybrid/Gaff.txt G.nobilis ~/Workspace/Robert/Gambusia/analysis/hybrid/West_Texas_nobilis.txt
cd ..
done
}
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "AxN_WTX.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "AxN_WTX.xlsx", overwrite = TRUE)
}
}

## G. affinis vs G. geiseri ##
{```{bash}```
#Generating data
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid
#Preparing dir
mkdir -p affinis/Ggei
cd affinis/Ggei
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
echo "Processing run $i"
cd run_$i
Rscript ~/bin/NewHybrids_vcf_to_txt.r /home/afields/Workspace/Robert/Gambusia/analysis/hybrid/Ga_Gg.recode.vcf newhybrids_input_AxG.txt
cd ..
done

#Adding header
for i in $(seq 1 5); do 
INDV=$(vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Ga_Gg.recode.vcf | wc -l);
sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_AxG.txt;
done

#Running iterations
for i in $(seq 1 5); do
cd run_$i
RAND1=$(echo $((RANDOM)))
RAND2=$(echo $((RANDOM)))
echo $RAND1 $RAND2 > starting_seeds.txt
newhybrids -d newhybrids_input_AxG.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui &
sleep 5
cd ..
done
}
#Formatting the output data
{```{bash}```
for i in $(seq 1 5); do
cd run_$i
vcfsamplenames ~/Workspace/Robert/Gambusia/analysis/hybrid/Ga_Gg.recode.vcf > tmp.indv
sed -i "1 s/^/Indv\n/" tmp.indv
paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
rm tmp.indv
Rscript ~/bin/NewHybrids_results_format.r G.affinis ~/Workspace/Robert/Gambusia/analysis/hybrid/Gaff.txt G.geiseri ~/Workspace/Robert/Gambusia/analysis/hybrid/Ggei.txt
cd ..
done
}
{```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "AxG.xlsx", colNames = T, sheetName="run1")
dat2 <- read.table("run_2/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat2$Max <- apply(dat2[,2:ncol(dat2)],1,max)
addWorksheet(wb, sheetName = "run2")
writeData(wb, "run2", dat2)
dat3 <- read.table("run_3/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat3$Max <- apply(dat3[,2:ncol(dat3)],1,max)
addWorksheet(wb, sheetName = "run3")
writeData(wb, "run3", dat3)
dat4 <- read.table("run_4/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat4$Max <- apply(dat4[,2:ncol(dat4)],1,max)
addWorksheet(wb, sheetName = "run4")
writeData(wb, "run4", dat4)
dat5 <- read.table("run_5/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat5$Max <- apply(dat5[,2:ncol(dat5)],1,max)
addWorksheet(wb, sheetName = "run5")
writeData(wb, "run5", dat5)
addWorksheet(wb, sheetName = "Summary")
saveWorkbook(wb, "AxG.xlsx", overwrite = TRUE)
}
}


### Adegenet simulations ###
## G. geseri vs New Mexico G. nobilis ##
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid
echo -e "Sample_ID\tSpecies"> species.txt
awk -v OFS="\t" '{print $0, "G.affinis"}' Gaff.txt >> species.txt
awk -v OFS="\t" '{print $0, "G.geiseri"}' Ggei.txt >> species.txt
cat Diamond_Y_nobilis.txt New_Mexico_nobilis.txt West_Texas_nobilis.txt | awk -v OFS="\t" '{print $0, "G.nobilis"}' >> species.txt
mkdir sims
cd sims
cp -s ../species.txt .
}
{```{R}```
#Loading Libraries
library(adegenet)
library(vcfR)
library(reshape2)
library(ggplot2)

#Importing data
vcf <- read.vcfR(file = "../Gg_NM.recode.vcf")
gen.vcf <- vcfR2genind(vcf)
strata <- read.table(file = "species.txt", header = TRUE, sep="\t")
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$Sample_ID),]
head(strata(gen.vcf))
rm(vcf)

#Define Strata
setPop(gen.vcf) <- ~Species
temp <- seppop(gen.vcf)
geis <- temp$`G.geiseri`
nob <- temp$`G.nobilis`

#Test for missing loci
geis.miss <- apply(geis@tab, 2, function(x) sum(is.na(x))/nInd(geis))
nob.miss <- apply(nob@tab, 2, function(x) sum(is.na(x))/nInd(nob))

rm.loci <- NULL
if(length(which(geis.miss  == 1))){rm.loci <- c(rm.loci, names(which(geis.miss == 1)))}
if(length(which(nob.miss  == 1))){rm.loci <- c(rm.loci, names(which(nob.miss == 1)))}
rm.loci <- unique(matrix(unlist(strsplit(rm.loci,"[.]")), ncol=2, byrow=T)[,1])

set.loc <- locNames(geis)[!locNames(geis) %in% rm.loci]
gen.vcf <- gen.vcf[, loc=set.loc]
geis <- geis[, loc=set.loc]
nob <- nob[, loc=set.loc]

#Simulate samples
F1 <- hybridize(geis, nob, n = 100, pop = "geisxnob")
geis_bx <- hybridize(geis, F1, n = 100, pop = "geis_bx")
nob_bx <- hybridize(nob, F1, n = 100, pop = "nob_bx")
F2 <- hybridize(F1, F1, n= 100, pop = "F1XF1")
pooled_gens <- repool(gen.vcf, F1, geis_bx, nob_bx, F2)

save(pooled_gens, file="G_x_N_NM_sims.gz", compress=T)
#load("G_x_N_NM_sims.gz")
#PCA
x <- scaleGen(pooled_gens, NA.method = "mean")

#pca with all loci
pca <- dudi.pca(x,center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

#plot2
col = c("#68a691", "#ddd392", "#ed6a5a", "#e59a2e", "#bdbdbd", "#636363", '#000000', '#799c41')
s.class(pca$li, pop(pooled_gens), xax=1, yax=2, cellipse=3, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0)
MIN.x <- par("usr")[1]*1.01
MAX.x <- par("usr")[2]
MIN.y <- par("usr")[3]
MAX.y <- par("usr")[4]

png("GgxGn_NM_sim.png", res=200, width=2000, height=2000)
s.class(pca$li[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1"), 1:2], pop(pooled_gens)[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1")], 
xax=1, yax=2, cellipse=10, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0, xlim=c(MIN.x,MAX.x), ylim=c(MIN.y,MAX.y))
points(pca$li[1:length(indNames(gen.vcf)),1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[1:length(indNames(gen.vcf))])])
legend(130, -100, legend=c("G. geiseri","G. nobilis"), col=col[1:2], pch=19, bty="n", cex=1, title=as.expression(bquote(bold("Species"))))
legend(130, -140, legend=c("F1 cross","G. geiseri BX","G. nobilis BX","F1 x F1"), col=col[3:7], pch=21, bty="n", cex=1, title = as.expression(bquote(bold("Simulation"))))
mtext("G. geiseri and G. nobilis (New Mexico)", 3, adj=0.05, cex=1.5)
text(-150, 30, "San Mateo River")
text(-155, -10, "West Texas")
#G.geiseri BX
points(pca$li["Gag-ES_MGL-16327", 1:2], pch=19, col=col[4])
points(pca$li["Gag-ES_MGL-16330", 1:2], pch=19, col=col[4])
#F1_cross
points(pca$li["Gnob_PSGn_2", 1:2], pch=19, col=col[3])
dev.off()
}}} #notepad cleanup

## G. geseri vs Diamond Y G. nobilis ##
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid/sims
}
{```{R}```
#Loading Libraries
library(adegenet)
library(vcfR)
library(reshape2)
library(ggplot2)

#Importing data
vcf <- read.vcfR(file = "../Gg_DY.recode.vcf")
gen.vcf <- vcfR2genind(vcf)
strata <- read.table(file = "species.txt", header = TRUE, sep="\t")
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$Sample_ID),]
head(strata(gen.vcf))
rm(vcf)

#Define Strata
setPop(gen.vcf) <- ~Species
temp <- seppop(gen.vcf)
geis <- temp$`G.geiseri`
nob <- temp$`G.nobilis`

#Test for missing loci
geis.miss <- apply(geis@tab, 2, function(x) sum(is.na(x))/nInd(geis))
nob.miss <- apply(nob@tab, 2, function(x) sum(is.na(x))/nInd(nob))

rm.loci <- NULL
if(length(which(geis.miss  == 1))){rm.loci <- c(rm.loci, names(which(geis.miss == 1)))}
if(length(which(nob.miss  == 1))){rm.loci <- c(rm.loci, names(which(nob.miss == 1)))}

if(length(rm.loci) > 0){
rm.loci <- unique(matrix(unlist(strsplit(rm.loci,"[.]")), ncol=2, byrow=T)[,1])
set.loc <- locNames(geis)[!locNames(geis) %in% rm.loci]
gen.vcf <- gen.vcf[, loc=set.loc]
geis <- geis[, loc=set.loc]
nob <- nob[, loc=set.loc]
}

#Simulate samples
F1 <- hybridize(geis, nob, n = 100, pop = "geisxnob")
geis_bx <- hybridize(geis, F1, n = 100, pop = "geis_bx")
nob_bx <- hybridize(nob, F1, n = 100, pop = "nob_bx")
F2 <- hybridize(F1, F1, n= 100, pop = "F1XF1")
pooled_gens <- repool(gen.vcf, F1, geis_bx, nob_bx, F2)

save(pooled_gens, file="G_x_N_DY_sims.gz", compress=T)
#load("G_x_N_DY_sims.gz")
#PCA
x <- scaleGen(pooled_gens, NA.method = "mean")

#pca with all loci
pca <- dudi.pca(x,center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

#plot2
col = c("#68a691", "#ddd392", "#ed6a5a", "#e59a2e", "#bdbdbd", "#636363", '#000000', '#799c41')
s.class(pca$li, pop(pooled_gens), xax=1, yax=2, cellipse=3, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0)
MIN.x <- par("usr")[1]
MAX.x <- par("usr")[2]
MIN.y <- par("usr")[3]
MAX.y <- par("usr")[4]

png("GgxGn_DY_sim.png", res=200, width=2000, height=2000)
s.class(pca$li[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1"), 1:2], pop(pooled_gens)[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1")], 
xax=1, yax=2, cellipse=10, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0, xlim=c(MIN.x,MAX.x), ylim=c(MIN.y,MAX.y))
points(pca$li[1:length(indNames(gen.vcf)),1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[1:length(indNames(gen.vcf))])])
legend(150, 150, legend=c("G. geiseri","G. nobilis"), col=col[1:2], pch=19, bty="n", cex=1, title=as.expression(bquote(bold("Species"))))
legend(150, 120, legend=c("F1 cross","G. geiseri BX","G. nobilis BX","F1 x F1"), col=col[3:7], pch=21, bty="n", cex=1, title = as.expression(bquote(bold("Simulation"))))
mtext("G. geiseri and G. nobilis (Diamond Y)", 3, adj=0.05, cex=1.5)
text(-196, -33, "San Mateo River")
text(-196, 20, "West Texas")
#G.geiseri BX
points(pca$li["Gag-ES_MGL-16327", 1:2], pch=19, col=col[4])
points(pca$li["Gag-ES_MGL-16330", 1:2], pch=19, col=col[4])
#F1_cross
points(pca$li["Gnob_PSGn_2", 1:2], pch=19, col=col[3])
dev.off()
}}}} #notepad cleanup

## G. geseri vs West Texas G. nobilis ##
{```{bash}```
cd /home/afields/Workspace/Robert/Gambusia/analysis/hybrid/sims
}
{```{R}```
#Loading Libraries
library(adegenet)
library(vcfR)
library(reshape2)
library(ggplot2)

#Importing data
vcf <- read.vcfR(file = "../Gg_WTX.recode.vcf")
gen.vcf <- vcfR2genind(vcf)
strata <- read.table(file = "species.txt", header = TRUE, sep="\t")
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$Sample_ID),]
head(strata(gen.vcf))
rm(vcf)

#Define Strata
setPop(gen.vcf) <- ~Species
temp <- seppop(gen.vcf)
geis <- temp$`G.geiseri`
nob <- temp$`G.nobilis`

#Test for missing loci
geis.miss <- apply(geis@tab, 2, function(x) sum(is.na(x))/nInd(geis))
nob.miss <- apply(nob@tab, 2, function(x) sum(is.na(x))/nInd(nob))

rm.loci <- NULL
if(length(which(geis.miss  == 1))){rm.loci <- c(rm.loci, names(which(geis.miss == 1)))}
if(length(which(nob.miss  == 1))){rm.loci <- c(rm.loci, names(which(nob.miss == 1)))}

if(length(rm.loci) > 0){
rm.loci <- unique(matrix(unlist(strsplit(rm.loci,"[.]")), ncol=2, byrow=T)[,1])
set.loc <- locNames(geis)[!locNames(geis) %in% rm.loci]
gen.vcf <- gen.vcf[, loc=set.loc]
geis <- geis[, loc=set.loc]
nob <- nob[, loc=set.loc]
}

#Simulate samples
F1 <- hybridize(geis, nob, n = 100, pop = "geisxnob")
geis_bx <- hybridize(geis, F1, n = 100, pop = "geis_bx")
nob_bx <- hybridize(nob, F1, n = 100, pop = "nob_bx")
F2 <- hybridize(F1, F1, n= 100, pop = "F1XF1")
pooled_gens <- repool(gen.vcf, F1, geis_bx, nob_bx, F2)

save(pooled_gens, file="G_x_N_WTX_sims.gz", compress=T)
#load("G_x_N_WTX_sims.gz")
#PCA
x <- scaleGen(pooled_gens, NA.method = "mean")

#pca with all loci
pca <- dudi.pca(x,center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

#plot2
col = c("#68a691", "#ddd392", "#ed6a5a", "#e59a2e", "#bdbdbd", "#636363", '#000000', '#799c41')
s.class(pca$li, pop(pooled_gens), xax=1, yax=2, cellipse=3, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0)
MIN.x <- par("usr")[1]
MAX.x <- par("usr")[2]
MIN.y <- par("usr")[3]
MAX.y <- par("usr")[4]

png("GgxGn_WTX_sim.png", res=200, width=2000, height=2000)
s.class(pca$li[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1"), 1:2], pop(pooled_gens)[pop(pooled_gens) %in% c("geisxnob", "geis_bx", "nob_bx", "F1XF1")], 
xax=1, yax=2, cellipse=10, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0, xlim=c(MIN.x,MAX.x), ylim=c(MIN.y,MAX.y))
points(pca$li[1:length(indNames(gen.vcf)),1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[1:length(indNames(gen.vcf))])])
legend(-150, -75, legend=c("G. geiseri","G. nobilis"), col=col[1:2], pch=19, bty="n", cex=1, title=as.expression(bquote(bold("Species"))))
legend(-150, -115, legend=c("F1 cross","G. geiseri BX","G. nobilis BX","F1 x F1"), col=col[3:7], pch=21, bty="n", cex=1, title = as.expression(bquote(bold("Simulation"))))
mtext("G. geiseri and G. nobilis (West Texas)", 3, adj=0.05, cex=1.5)
text(160, 70,"Phantom Lake")
text(160, -35,"Balmorhea and East Sadia")
#G.geiseri BX
points(pca$li["Gag-ES_MGL-16327", 1:2], pch=19, col=col[4])
points(pca$li["Gag-ES_MGL-16330", 1:2], pch=19, col=col[4])
#F1_cross
points(pca$li["Gnob_PSGn_2", 1:2], pch=19, col=col[3])
#G.nobilis BX
points(pca$li["Gnob_PSGn_5", 1:2], pch=19, col=col[5])
points(pca$li["Gnob_PSGn_8", 1:2], pch=19, col=col[5])
dev.off()
}}}} #notepad cleanup


