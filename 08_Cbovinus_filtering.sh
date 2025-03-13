#Filtering the Cyprinodon bovinus TotalRawSNPs.vcf
#Charting the raw data
```{bash}```
sed -i 's/-RG//g' popmap
charts.sh TotalRawSNPs.vcf &

#Split Multiallelic Freebayes calls
```{bash}```
vcfallelicprimitives -k -g TotalRawSNPs.vcf | sed 's:\.|\.:\.\/\.:g' > TRS.prim

#Remove indels
```{bash}```
vcftools --vcf TRS.prim --recode-INFO-all --recode --out SNP.TRS --remove-indels

#Basic quality filters of data (Minimum of 2 alleles, Depth of at least 5, Minimum quality of at least 20)
```{bash}```
vcftools --vcf SNP.TRS.recode.vcf --out SNP.TRS.QC --recode --recode-INFO-all --min-alleles 2 --minDP 5 --minQ 20 --min-meanDP 20

#Filters allelic balance, quality vs depth, strand representation and paired read representation
```{bash}```
dDocent_filters SNP.TRS.QC.recode.vcf SNP.TRS.dDocent
## Input
#no
#100000

#Filtering singleton and doubleton loci for depth (Singletons >20 reads; Doubletons > 10 reads)
```{bash}```
vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out out --singletons

awk ' $3=="S" {print $1, $2}' out.singletons > sing.loci
awk ' $3=="D" {print $1, $2}' out.singletons > doub.loci

vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out SNP.TRS.F01a --recode --recode-INFO-all --exclude-positions sing.loci
vcftools --vcf SNP.TRS.F01a.recode.vcf --out SNP.TRS.F01b --recode --recode-INFO-all --exclude-positions doub.loci

vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out SNP.TRS.F01.sing --recode --recode-INFO-all --positions sing.loci
vcftools --vcf SNP.TRS.F01.sing.recode.vcf --out SNP.TRS.F02.sing --recode --recode-INFO-all --minDP 20

vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out SNP.TRS.F01.doub --recode --recode-INFO-all --positions doub.loci
vcftools --vcf SNP.TRS.F01.doub.recode.vcf --out SNP.TRS.F02.doub --recode --recode-INFO-all --minDP 10

vcf-concat SNP.TRS.F02.sing.recode.vcf SNP.TRS.F02.doub.recode.vcf SNP.TRS.F01b.recode.vcf > SNP.TRS.F02.vcf

rm out.singletons

#Filter loci that have high variation in depth across a locus with an individual
```{bash}```
vcftools --vcf SNP.TRS.F02.vcf --out out --geno-depth

```{R}```
gdepth<-read.table(file="out.gdepth", head=T)
gdepth[gdepth==-1]<-NA

for (i in 3:dim(gdepth)[2]) {
temp<-aggregate(gdepth[,i],by=list(gdepth[,1]), sd)
if(i==3){indv.site.sd<-data.frame(temp,row.names=1)} else
{indv.site.sd[,(i-2)]<-temp[,2]}
}
colnames(indv.site.sd)<-colnames(gdepth[3:dim(gdepth)[2]])
tmp<-apply(indv.site.sd, 1, mean, na.rm=T)
tmp2<-unique(c(names(which(tmp>25))))
length(tmp)
length(tmp2)
write.table(tmp2,file="bad.loci.sd", quote=F, col.names=F, row.names=F)
q("no")

```{bash}```
grep "dDocent" SNP.TRS.F02.vcf | cut -f 1,2 | uniq | tail -n +2 > contigs.txt
grep -wf bad.loci.sd contigs.txt > bad.loci
vcftools --vcf SNP.TRS.F02.vcf  --out SNP.TRS.F03 --exclude-positions bad.loci --recode-INFO-all --recode

charts.sh SNP.TRS.F03.recode.vcf &

#Filter data using decision tree
```{bash}```
mkdir SOL
cp SNP.TRS.F03.recode.vcf SOL/
cd SOL
SOL.filter1.sh SNP.TRS.F03.recode.vcf

cp vcf/A.1.1.1.2.2.3.1.2.SNP.finalb.1.recode.vcf ../
cd ..
charts.sh A.1.1.1.2.2.3.1.2.SNP.finalb.1.recode.vcf

#Removing loci with super high depth
```{bash}```
dDocent_filters A.1.1.1.2.2.3.1.2.SNP.finalb.1.recode.vcf SNP.TRS.F04
##Input
#no
#no

charts.sh SNP.TRS.F04.FIL.recode.vcf

#Removing paralogs and indivduals with high missing data#
```{bash}```
mkdir haplotyping
cd haplotyping

cp -s ../../reference/Cbovinus_reference.fasta ./reference.fasta
cp -s ../popmap .
cp -s ../mapping/*-RG.bam .
rm cat-RRG.bam*
ls *.bam | while read i; do echo "Processing $i"; samtools index $i; done
cp -s ../SNP.TRS.F04.FIL.recode.vcf .

rad_haplotyper.pl -v SNP.TRS.F04.FIL.recode.vcf -r reference.fasta -p popmap -x 20 -m 0.8 -o SNP.TRS.F05.vcf -g SNP.TRS.F05.gen

cd ..
cp haplotyping/SNP.TRS.F05.vcf .
charts.sh SNP.TRS.F05.vcf
