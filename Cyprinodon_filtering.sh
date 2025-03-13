#Filtering the Cyprinodon TotalRawSNPs.vcf
#Charting the raw data
{```{bash}`
sed -i 's/-RG//g' popmap
charts.sh TotalRawSNPs.vcf &
```}

#Split Multiallelic Freebayes calls
{```{bash}
vcfallelicprimitives -k -g TotalRawSNPs.vcf | sed 's:\.|\.:\.\/\.:g' > TRS.prim
```}

#Remove indels
{```{bash}
vcftools --vcf TRS.prim --recode-INFO-all --recode --out SNP.TRS --remove-indels
```}
#After filtering, kept 118 out of 118 Individuals
#After filtering, kept 1706116 out of a possible 1834508 Sites
#Run Time = 1125.00 seconds

#Basic quality filters of data (Minimum of 2 alleles, Depth of at least 5, Minimum quality of at least 20)
{```{bash}
vcftools --vcf SNP.TRS.recode.vcf --out SNP.TRS.QC --recode --recode-INFO-all --min-alleles 2 --minDP 5 --minQ 20 --min-meanDP 20
```}
#After filtering, kept 118 out of 118 Individuals
#After filtering, kept 258397 out of a possible 1706116 Sites
#Run Time = 151.00 seconds

#Filters allelic balance, quality vs depth, strand representation and paired read representation
{```{bash}
dDocent_filters SNP.TRS.QC.recode.vcf SNP.TRS.dDocent
#no
#100000
```}
#Remaining sites
# 164605

#Filtering singleton and doubleton loci for depth (Singletons >20 reads; Doubletons > 10 reads)
{```{bash}
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
```}

#Filter loci that have high variation in depth across a locus with an individual
```{bash}
vcftools --vcf SNP.TRS.F02.vcf --out out --geno-depth

R
```

```{R}
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
```

```{bash}
grep "dDocent" SNP.TRS.F02.vcf | cut -f 1,2 | uniq | tail -n +2 > contigs.txt
grep -wf bad.loci.sd contigs.txt > bad.loci
vcftools --vcf SNP.TRS.F02.vcf  --out SNP.TRS.F03 --exclude-positions bad.loci --recode-INFO-all --recode
#After filtering, kept 118 out of 118 Individuals
#After filtering, kept 163267 out of a possible 164624 Sites
#Run Time = 94.00 seconds

charts.sh SNP.TRS.F03.recode.vcf &
```

#Filter data using decision tree
```{bash}
mkdir SOL
cp SNP.TRS.F03.recode.vcf SOL/
cd SOL
SOL.filter1.sh SNP.TRS.F03.recode.vcf

#The Final_Filter.count in the SOL folder needs to be exported and evaluated to determine the best filtering path
#Commonly B.1.3.3.2.3.2.1.2.SNP.finalb.1.recode.vcf has the best compromise keeping the most loci and the most individuals

cp vcf/A.1.1.1.2.2.3.1.2.SNP.finalb.1.recode.vcf ../
cd ..
charts.sh A.1.1.1.2.2.3.1.2.SNP.finalb.1.recode.vcf
```

#Removing loci with super high depth
```{bash}
dDocent_filters A.1.1.1.2.2.3.1.2.SNP.finalb.1.recode.vcf SNP.TRS.F04
#no
#no

charts.sh SNP.TRS.F04.FIL.recode.vcf
```
#Remaining sites
# 33021

#Filtering individuals with > 20% missing data
```{bash}
vcftools --vcf SNP.TRS.F04.FIL.recode.vcf --out SNP.TRS.F04.FIL --missing-indv
mawk -v x=0.2 '$5 > x' SNP.TRS.F04.FIL.imiss | cut -f1 > lowDP.indv
vcftools --vcf SNP.TRS.F04.FIL.recode.vcf --out SNP.TRS.F05 --remove lowDP.indv --recode --recode-INFO-all
#After filtering, kept 113 out of 115 Individuals
#After filtering, kept 33021 out of a possible 33021 Sites
#Run Time = 36.00 seconds

charts.sh SNP.TRS.F05.recode.vcf
```

#Removing paralogs and indivduals with high missing data#
```{bash}
#Earth
mkdir haplotyping
cd haplotyping

#HPC
scp reference/Cyprinodon_reference.fasta afields@earth.ad.tamucc.edu:/home/afields/Workspace/Robert/Cyprinodon/reference/
scp mapping/I*/*.bam afields@earth.ad.tamucc.edu:/home/afields/Workspace/Robert/Cyprinodon/filtering/haplotyping

#Earth
cp -s ../../reference/Cyprinodon_reference.fasta ./reference.fasta
cp -s ../popmap .
rm cat-RRG.bam*
ls *.bam | while read i; do echo "Processing $i"; samtools index $i; done
cp -s ../SNP.TRS.F05.recode.vcf .

rad_haplotyper.pl -v SNP.TRS.F05.recode.vcf -r reference.fasta -p popmap -x 20 -m 0.8 -o SNP.TRS.F06.vcf -g SNP.TRS.F06.gen
````

#Preparing analysis directory
{```{bash}```
mkdir ../analysis
cd ../analysis
cp ../filtering/SNP.TRS.F06.recode.vcf ../filtering/haplotyping/SNP.TRS.F06.gen .
}

#strata updated by hand in excel
#indv.csv in SNP.TRS.F04.FIL.2023-01-09 updated
}