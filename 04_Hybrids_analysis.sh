### Hybrid analysis of the Cyprinodon ###

## NewHybrids ##
#Taking files with individuals listed by group and combining them for pairwise analysis
```{bash}```
#Make bovinus and variegatus lists
vcfsamplenames SNP.TRS.F06.vcf > Cbov_Cvar.txt
grep Cyb Cbov_Cvar.txt > Cbov.txt
grep Cyv Cbov_Cvar.txt > Cvar.txt

#Generating data
```{bash}```
#Preparing dir
mkdir new_hybrids && cd new_hybrids
for i in $(seq 1 5); do mkdir run_$i; done

for i in $(seq 1 5); do 
   echo "Processing run $i"
   cd run_$i
   Rscript ~/bin/NewHybrids_vcf_to_txt.r ../../SNP.TRS.F06.vcf newhybrids_input_bovxvar.txt
   cd ..
done

#Adding header
```{bash}```
for i in $(seq 1 5); do 
   INDV=$(vcfsamplenames SNP.TRS.F06.vcf | wc -l);
   sed -i "1 s/^/NumIndivs $INDV\nNumLoci 150\nDigits 1\nFormat Lumped\n\n/" run_$i/newhybrids_input_bovxvar.txt;
done

#Running iterations
```{bash}```
for i in $(seq 1 5); do
   cd run_$i
   RAND1=$(echo $((RANDOM)))
   RAND2=$(echo $((RANDOM)))
   echo $RAND1 $RAND2 > starting_seeds.txt
   newhybrids -d newhybrids_input_bovxvar.txt -c ~/bin/genotype_frequencies.txt -s $RAND1 $RAND2 --no-gui &
   sleep 5
   cd ..
done
}

#Formatting the output data
```{bash}```
for i in $(seq 1 5); do
   cd run_$i
   vcfsamplenames ../SNP.TRS.F06.vcf > tmp.indv
   sed -i "1 s/^/Indv\n/" tmp.indv
   paste tmp.indv aa-ScaledLikelihood.txt > aa-ScaledLikelihood_indv.txt
   rm tmp.indv
   Rscript ~/bin/NewHybrids_results_format.r C.bovinus ../Cbov.txt C.variegatus ../Cvar.txt
   cd ..
done

```{R}```
library(openxlsx)
dat <- read.table("run_1/aa-ScaledLikelihood_edit.txt", head=T, sep="\t")
dat$Max <- apply(dat[,2:ncol(dat)],1,max)
wb <- write.xlsx(dat, "bov_var.xlsx", colNames = T, sheetName="run1")
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
saveWorkbook(wb, "bov_var.xlsx", overwrite = TRUE)

### Adegenet simulations ###
```{bash}```
echo -e "Sample_ID\tSpecies"> species.txt
awk -v OFS="\t" '{print $0, "C.bovinus"}' Cbov.txt >> species.txt
awk -v OFS="\t" '{print $0, "C.variegatus"}' Cvar.txt >> species.txt

mkdir sims
cd sims
cp -s ../species.txt .

```{R}```
#Loading Libraries
library(adegenet)
library(vcfR)
library(reshape2)
library(ggplot2)

#Importing data
```{R}```
vcf <- read.vcfR(file = "../SNP.TRS.F06.vcf")
gen.vcf <- vcfR2genind(vcf)
strata <- read.table(file = "species.txt", header = TRUE, sep="\t")
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$Sample_ID),]
head(strata(gen.vcf))
rm(vcf)

#Define Strata
```{R}```
setPop(gen.vcf) <- ~Species
tmp <- seppop(gen.vcf)
bov <- tmp$`C.bovinus`
var <- tmp$`C.variegatus`

#Test for missing loci
```{R}```
bov.miss <- apply(bov@tab, 2, function(x) sum(is.na(x))/nInd(bov))
var.miss <- apply(var@tab, 2, function(x) sum(is.na(x))/nInd(var))

rm.loci <- NULL
if(length(which(bov.miss  == 1))){rm.loci <- c(rm.loci, names(which(bov.miss == 1)))}
if(length(which(var.miss  == 1))){rm.loci <- c(rm.loci, names(which(var.miss == 1)))}
rm.loci <- unique(matrix(unlist(strsplit(rm.loci,"[.]")), ncol=2, byrow=T)[,1])

set.loc <- locNames(bov)[!locNames(bov) %in% rm.loci]
gen.vcf <- gen.vcf[, loc=set.loc]
bov <- bov[, loc=set.loc]
var <- var[, loc=set.loc]

#Simulate samples
```{R}```
F1 <- hybridize(bov, var, n = 100, pop = "bovxvar")
bov_bx <- hybridize(bov, F1, n = 100, pop = "bov_bx")
var_bx <- hybridize(var, F1, n = 100, pop = "var_bx")
F2 <- hybridize(F1, F1, n= 100, pop = "F1XF1")
pooled_gens <- repool(gen.vcf, F1, bov_bx, var_bx, F2)

save(pooled_gens, file="bov_x_var_NM_sims.gz", compress=T)

#PCA
```{R}```
x <- scaleGen(pooled_gens, NA.method = "mean")
pca <- dudi.pca(x,center = FALSE, scale = FALSE, scannf = FALSE, nf = 4)

#Plot parameters
```{R}```
col = c("#68a691", "#ddd392", "#ed6a5a", "#e59a2e", "#bdbdbd", "#636363", '#000000', '#799c41')
s.class(pca$li, pop(pooled_gens), xax=1, yax=2, cellipse=3, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0)
MIN.x <- par("usr")[1]
MAX.x <- par("usr")[2]
MIN.y <- par("usr")[3]
MAX.y <- par("usr")[4]

#plot
```{R}```
png("Cbov_x_Cvar_sim.png", res=200, width=2000, height=2000)
s.class(pca$li[pop(pooled_gens) %in% c("bovxvar", "bov_bx", "var_bx", "F1XF1"), 1:2], pop(pooled_gens)[pop(pooled_gens) %in% c("bovxvar", "bov_bx", "var_bx", "F1XF1")], 
xax=1, yax=2, cellipse=10, col=transp(col, 1), origin=c(0,0), axesell = FALSE, cstar=0, cpoint=0, grid=FALSE, addaxes=TRUE, clabel=0, xlim=c(MIN.x,MAX.x), ylim=c(MIN.y,MAX.y))
points(pca$li[1:length(indNames(gen.vcf)),1:2], pch=19, col=col[as.numeric(pop(pooled_gens)[1:length(indNames(gen.vcf))])])
legend(130, -100, legend=c("C. bovinus","C. variegatus"), col=col[1:2], pch=19, bty="n", cex=1, title=as.expression(bquote(bold("Species"))))
legend(130, -140, legend=c("F1 cross","C. bovinus BX","C. variegatus BX","F1 x F1"), col=col[3:7], pch=21, bty="n", cex=1, title = as.expression(bquote(bold("Simulation"))))
mtext("C. bovinus and C. variegatus", 3, adj=0.05, cex=1.5)
dev.off()
