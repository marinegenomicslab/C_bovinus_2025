#!/bin/bash

# dDocent filters (requires user input) - will stop script
dDocent_filters B.1.3.3.2.3.2.1.2.recode.vcf B.1.3.3.2.3.2.1.2 PE no
dDocent_filters B.1.3.3.2.3.2.2.1.recode.vcf B.1.3.3.2.3.2.2.2 PE no
dDocent_filters A.1.1.1.2.2.3.1.2.recode.vcf A.1.1.1.2.2.3.1.2 PE no
dDocent_filters B.1.3.3.2.3.4.1.2.recode.vcf B.1.3.3.2.3.4.1.2 PE no
dDocent_filters B.1.3.3.2.3.4.2.1.recode.vcf B.1.3.3.2.3.4.2.1 PE no
# Remove INDELS
vcfallelicprimitives B.1.3.3.2.3.2.1.2.recode.vcf --keep-info --keep-geno > prim.vcf 
vcftools --vcf prim.vcf --out B.1.3.3.2.3.2.1.2.SNP --remove-indels --recode --recode-INFO-all 
vcfallelicprimitives B.1.3.3.2.3.2.2.2.recode.vcf --keep-info --keep-geno > prim.vcf 
vcftools --vcf prim.vcf --out B.1.3.3.2.3.2.2.2.SNP --remove-indels --recode --recode-INFO-all 
vcfallelicprimitives A.1.1.1.2.2.3.1.2.recode.vcf --keep-info --keep-geno > prim.vcf 
vcftools --vcf prim.vcf --out A.1.1.1.2.2.3.1.2.SNP --remove-indels --recode --recode-INFO-all 
vcfallelicprimitives B.1.3.3.2.3.4.1.2.recode.vcf --keep-info --keep-geno > prim.vcf 
vcftools --vcf prim.vcf --out B.1.3.3.2.3.4.1.2.SNP --remove-indels --recode --recode-INFO-all 
vcfallelicprimitives B.1.3.3.2.3.4.2.1.recode.vcf --keep-info --keep-geno > prim.vcf 
vcftools --vcf prim.vcf --out B.1.3.3.2.3.4.2.1.SNP --remove-indels --recode --recode-INFO-all
# Final Threshold values Filter MAF
vcftools --vcf B.1.3.3.2.3.2.1.2.SNP.recode.vcf --out B.1.3.3.2.3.2.1.2.SNP.final0 --maf 0.05 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.2.2.2.SNP.recode.vcf --out B.1.3.3.2.3.2.2.2.SNP.final0 --maf 0.05 --recode --recode-INFO-all 
vcftools --vcf A.1.1.1.2.2.3.1.2.SNP.recode.vcf --out A.1.1.1.2.2.3.1.2.SNP.final0 --maf 0.05 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.4.1.2.SNP.recode.vcf --out B.1.3.3.2.3.4.1.2.SNP.final0 --maf 0.05 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.4.2.1.SNP.recode.vcf --out B.1.3.3.2.3.4.2.1.SNP.final0 --maf 0.05 --recode --recode-INFO-all
# filter missing loci
vcftools --vcf B.1.3.3.2.3.2.1.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.2.1.2.finala.1 --max-missing 0.95 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.2.1.2.finala.1.recode.vcf --out B_1.3.3.2.3.2.1.2.finala.1 --missing-indv 
vcftools --vcf B.1.3.3.2.3.2.2.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.2.2.2.finala.1 --max-missing 0.95 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.2.2.2.finala.1.recode.vcf --out B_1.3.3.2.3.2.2.2.finala.1 --missing-indv 
vcftools --vcf A.1.1.1.2.2.3.1.2.SNP.final0.recode.vcf --out A.1.1.1.2.2.3.1.2.finala.1 --max-missing 0.95 --recode --recode-INFO-all 
vcftools --vcf A.1.1.1.2.2.3.1.2.finala.1.recode.vcf --out A_1.1.1.2.2.3.1.2.finala.1 --missing-indv 
vcftools --vcf B.1.3.3.2.3.4.1.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.4.1.2.finala.1 --max-missing 0.95 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.4.1.2.finala.1.recode.vcf --out B_1.3.3.2.3.4.1.2.finala.1 --missing-indv 
vcftools --vcf B.1.3.3.2.3.4.2.1.SNP.final0.recode.vcf --out B.1.3.3.2.3.4.2.1.finala.1 --max-missing 0.95 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.4.2.1.finala.1.recode.vcf --out B_1.3.3.2.3.4.2.1.finala.1 --missing-indv 
vcftools --vcf B.1.3.3.2.3.2.1.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.2.1.2.finala.2 --max-missing 0.9 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.2.1.2.finala.2.recode.vcf --out B_1.3.3.2.3.2.1.2.finala.2 --missing-indv 
vcftools --vcf B.1.3.3.2.3.2.2.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.2.2.2.finala.2 --max-missing 0.9 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.2.2.2.finala.2.recode.vcf --out B_1.3.3.2.3.2.2.2.finala.2 --missing-indv 
vcftools --vcf A.1.1.1.2.2.3.1.2.SNP.final0.recode.vcf --out A.1.1.1.2.2.3.1.2.finala.2 --max-missing 0.9 --recode --recode-INFO-all 
vcftools --vcf A.1.1.1.2.2.3.1.2.finala.2.recode.vcf --out A_1.1.1.2.2.3.1.2.finala.2 --missing-indv 
vcftools --vcf B.1.3.3.2.3.4.1.2.SNP.final0.recode.vcf --out B.1.3.3.2.3.4.1.2.finala.2 --max-missing 0.9 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.4.1.2.finala.2.recode.vcf --out B_1.3.3.2.3.4.1.2.finala.2 --missing-indv 
vcftools --vcf B.1.3.3.2.3.4.2.1.SNP.final0.recode.vcf --out B.1.3.3.2.3.4.2.1.finala.2 --max-missing 0.9 --recode --recode-INFO-all 
vcftools --vcf B.1.3.3.2.3.4.2.1.finala.2.recode.vcf --out B_1.3.3.2.3.4.2.1.finala.2 --missing-indv
# filter missing individuals
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.1.2.finala.1.imiss | cut -f1 > lowDP.indv 
vcftools --vcf B.1.3.3.2.3.2.1.2.finala.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.2.SNP.finalb.1 
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.2.2.finala.1.imiss | cut -f1 > lowDP.indv 
vcftools --vcf B.1.3.3.2.3.2.2.2.finala.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.2.SNP.finalb.1 
mawk -v x=0.5 '$5 > x' A_1.1.1.2.2.3.1.2.finala.1.imiss | cut -f1 > lowDP.indv 
vcftools --vcf A.1.1.1.2.2.3.1.2.finala.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.2.3.1.2.SNP.finalb.1 
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.4.1.2.finala.1.imiss | cut -f1 > lowDP.indv 
vcftools --vcf B.1.3.3.2.3.4.1.2.finala.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.1.2.SNP.finalb.1 
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.4.2.1.finala.1.imiss | cut -f1 > lowDP.indv 
vcftools --vcf B.1.3.3.2.3.4.2.1.finala.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.2.1.SNP.finalb.1 
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.1.2.finala.2.imiss | cut -f1 > lowDP.indv 
vcftools --vcf B.1.3.3.2.3.2.1.2.finala.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.2.SNP.finalb.2 
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.2.2.finala.2.imiss | cut -f1 > lowDP.indv 
vcftools --vcf B.1.3.3.2.3.2.2.2.finala.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.2.SNP.finalb.2 
mawk -v x=0.5 '$5 > x' A_1.1.1.2.2.3.1.2.finala.2.imiss | cut -f1 > lowDP.indv 
vcftools --vcf A.1.1.1.2.2.3.1.2.finala.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.2.3.1.2.SNP.finalb.2 
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.4.1.2.finala.2.imiss | cut -f1 > lowDP.indv 
vcftools --vcf B.1.3.3.2.3.4.1.2.finala.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.1.2.SNP.finalb.2 
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.4.2.1.finala.2.imiss | cut -f1 > lowDP.indv 
vcftools --vcf B.1.3.3.2.3.4.2.1.finala.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.2.1.SNP.finalb.2
# remove unnecessary files
rm *.vcfidx *.count *.DEPTH *.depth meandepthpersite *.frq *.imiss lowDP.indv *.ldepth *.lmiss *.lowQDloci *.qual A_*.log B_*.log
# create folder with all VCF-files
mkdir vcf 
mv *.vcf vcf/
# create folder with all log-files
mkdir Filter.logs 
mv *.log Filter.logs/
# delete rmaining files
rm B.* A.*
# generate filter stats (SNP, Contig, Indv)
cd vcf 
echo "FILTER SNP CONTIG INDV" > Final_Filter.count 
for i in *.vcf 
do
  SNP=$(grep -cv '#' $i)
  CONTIG=$(grep -v '#' $i | cut -f 1 | sort | uniq | wc -l)
  INDV=$(vcfsamplenames $i | wc -l)
  echo "$i $SNP $CONTIG $INDV" >> ../Final_Filter.count 
done
# SOL
