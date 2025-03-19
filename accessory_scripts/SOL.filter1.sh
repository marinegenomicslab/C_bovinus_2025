#!/bin/bash

#Usage:
#SOL.filter1.sh <file to analyzed>

# FILTER 0: SEQUENCE QUALITY
vcftools --vcf $1 --out F0 --minQ 20 --recode --recode-INFO-all
vcftools --vcf $1 --out F0a --minQ 10 --recode --recode-INFO-all
vcftools --vcf $1 --out F0b --minQ 15 --recode --recode-INFO-all
# FILTER 1: GENOTYPE CALL RATE & MINIMUM ALLELE COUNT
vcftools --vcf F0.recode.vcf --out A --max-missing 0.3 --mac 3 --recode --recode-INFO-all
vcftools --vcf F0.recode.vcf --out B0 --max-missing 0.5 --mac 5 --recode --recode-INFO-all
vcftools --vcf F0.recode.vcf --out B --max-missing 0.5 --mac 3 --recode --recode-INFO-all
# FILTER 2: ELIMINATE SITES BY REQUIRED MINIMUM DEPTH
vcftools --vcf A.recode.vcf --out A.1 --minDP 3 --recode --recode-INFO-all
vcftools --vcf A.recode.vcf --out A.2 --minDP 5 --recode --recode-INFO-all
vcftools --vcf A.recode.vcf --out A.3 --minDP 10 --recode --recode-INFO-all
vcftools --vcf B.recode.vcf --out B.1 --minDP 3 --recode --recode-INFO-all
vcftools --vcf B.recode.vcf --out B.2 --minDP 5 --recode --recode-INFO-all
vcftools --vcf B.recode.vcf --out B.3 --minDP 10 --recode --recode-INFO-all
# FILTER 3: ELIMINATE INDIVIDUALS/LOCI WITH HIGH PROPORTION OF MISSING DATA I
# A: filter loci
vcftools --vcf A.1.recode.vcf --out A.1a --max-missing 0.5 --recode --recode-INFO-all
vcftools --vcf A.1a.recode.vcf --out A_1a --missing-indv
vcftools --vcf A.2.recode.vcf --out A.2a --max-missing 0.5 --recode --recode-INFO-all
vcftools --vcf A.2a.recode.vcf --out A_2a --missing-indv
vcftools --vcf A.3.recode.vcf --out A.3a --max-missing 0.5 --recode --recode-INFO-all
vcftools --vcf B.1.recode.vcf --out B.1a --max-missing 0.6 --recode --recode-INFO-all
vcftools --vcf B.1a.recode.vcf --out B_1a --missing-indv
vcftools --vcf B.2.recode.vcf --out B.2a --max-missing 0.6 --recode --recode-INFO-all
vcftools --vcf B.2a.recode.vcf --out B_2a --missing-indv
# B: filter individuals
mawk -v x=0.99 '$5 > x' A_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1
mawk -v x=0.95 '$5 > x' A_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.2
mawk -v x=0.9 '$5 > x' A_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.3
mawk -v x=0.9 '$5 > x' A_2a.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.2a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.2.1
mawk -v x=0.99 '$5 > x' B_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.1
mawk -v x=0.95 '$5 > x' B_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.2
mawk -v x=0.9 '$5 > x' B_1a.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1a.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3
# FILTER 4: MINOR ALLELE FREQUENCY AND MEAN DEPTH I
vcftools --vcf A.1.1.recode.vcf --out A.1.1.1 --maf 0.01 --min-meanDP 5 --recode --recode-INFO-all
vcftools --vcf A.1.1.recode.vcf --out A.1.1.2 --maf 0.02 --min-meanDP 5 --recode --recode-INFO-all
vcftools --vcf A.1.1.recode.vcf --out A.1.1.3 --maf 0.01 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf A.1.1.recode.vcf --out A.1.1.4 --maf 0.02 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf B.1.3.recode.vcf --out B.1.3.1 --maf 0.01 --min-meanDP 5 --recode --recode-INFO-all
vcftools --vcf B.1.3.recode.vcf --out B.1.3.2 --maf 0.01 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf B.1.3.recode.vcf --out B.1.3.3 --maf 0.02 --min-meanDP 5 --recode --recode-INFO-all
vcftools --vcf B.1.3.recode.vcf --out B.1.3.4 --maf 0.02 --min-meanDP 10 --recode --recode-INFO-all
# FILTER 5: ELIMINATE INDIVIDUALS/LOCI WITH HIGH PROPORTION OF MISSING DATA II
# A: filter loci
vcftools --vcf A.1.1.1.recode.vcf --out A.1.1.1.1 --max-missing 0.75 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.1.recode.vcf --out A_1.1.1.1 --missing-indv
vcftools --vcf A.1.1.1.recode.vcf --out A.1.1.1.2 --max-missing 0.8 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.recode.vcf --out A_1.1.1.2 --missing-indv
vcftools --vcf B.1.3.3.recode.vcf --out B.1.3.3.1 --max-missing 0.75 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.1.recode.vcf --out B_1.3.3.1 --missing-indv
vcftools --vcf B.1.3.3.recode.vcf --out B.1.3.3.2 --max-missing 0.8 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.recode.vcf --out B_1.3.3.2 --missing-indv
vcftools --vcf B.1.3.3.recode.vcf --out B.1.3.3.3 --max-missing 0.85 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.3.recode.vcf --out B_1.3.3.3 --missing-indv
# B: filter individuals
mawk -v x=0.95 '$5 > x' A_1.1.1.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.1.1
mawk -v x=0.9 '$5 > x' A_1.1.1.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.1
mawk -v x=0.95 '$5 > x' A_1.1.1.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.2
mawk -v x=0.85 '$5 > x' B_1.3.3.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.1
mawk -v x=0.8 '$5 > x' B_1.3.3.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.2
mawk -v x=0.75 '$5 > x' B_1.3.3.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3
mawk -v x=0.85 '$5 > x' B_1.3.3.3.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.3.1
mawk -v x=0.8 '$5 > x' B_1.3.3.3.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.3.2
mawk -v x=0.75 '$5 > x' B_1.3.3.3.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.3.3
# FILTER 6: MINOR ALLELE FREQUENCY AND MEAN DEPTH II
vcftools --vcf A.1.1.1.1.1.recode.vcf --out A.1.1.1.1.1.1 --maf 0.03 --min-meanDP 20 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.1.1.recode.vcf --out A.1.1.1.1.1.2 --maf 0.05 --min-meanDP 20 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.1 --maf 0.05 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.2 --maf 0.05 --min-meanDP 15 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.3 --maf 0.05 --min-meanDP 20 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.4 --maf 0.03 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.5 --maf 0.03 --min-meanDP 15 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.recode.vcf --out A.1.1.1.2.2.6 --maf 0.03 --min-meanDP 20 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.recode.vcf --out B.1.3.3.2.3.1 --maf 0.03 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.recode.vcf --out B.1.3.3.2.3.2 --maf 0.03 --min-meanDP 20 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.recode.vcf --out B.1.3.3.2.3.3 --maf 0.05 --min-meanDP 10 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.recode.vcf --out B.1.3.3.2.3.4 --maf 0.05 --min-meanDP 20 --recode --recode-INFO-all
# FILTER 7: ELIMINATE INDIVIDUALS/LOCI WITH HIGH PROPORTION OF MISSING DATA III
# A: filter loci
vcftools --vcf A.1.1.1.2.2.3.recode.vcf --out A.1.1.1.2.2.3.1 --max-missing 0.85 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.3.1.recode.vcf --out A_1.1.1.2.2.3.1 --missing-indv
vcftools --vcf A.1.1.1.2.2.3.recode.vcf --out A.1.1.1.2.2.3.2 --max-missing 0.9 --recode --recode-INFO-all
vcftools --vcf A.1.1.1.2.2.3.2.recode.vcf --out A_1.1.1.2.2.3.2 --missing-indv
vcftools --vcf B.1.3.3.2.3.2.recode.vcf --out B.1.3.3.2.3.2.1 --max-missing 0.85 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.2.1.recode.vcf --out B_1.3.3.2.3.2.1 --missing-indv
vcftools --vcf B.1.3.3.2.3.2.recode.vcf --out B.1.3.3.2.3.2.2 --max-missing 0.9 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.2.2.recode.vcf --out B_1.3.3.2.3.2.2 --missing-indv
vcftools --vcf B.1.3.3.2.3.4.recode.vcf --out B.1.3.3.2.3.4.1 --max-missing 0.85 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.4.1.recode.vcf --out B_1.3.3.2.3.4.1 --missing-indv
vcftools --vcf B.1.3.3.2.3.4.recode.vcf --out B.1.3.3.2.3.4.2 --max-missing 0.9 --recode --recode-INFO-all
vcftools --vcf B.1.3.3.2.3.4.2.recode.vcf --out B_1.3.3.2.3.4.2 --missing-indv
# B: filter individuals
mawk -v x=0.85 '$5 > x' A_1.1.1.2.2.3.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.2.2.3.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.2.3.1.1
mawk -v x=0.8 '$5 > x' A_1.1.1.2.2.3.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf A.1.1.1.2.2.3.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out A.1.1.1.2.2.3.1.2
mawk -v x=0.7 '$5 > x' B_1.3.3.2.3.2.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.1
mawk -v x=0.65 '$5 > x' B_1.3.3.2.3.2.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.2
mawk -v x=0.6 '$5 > x' B_1.3.3.2.3.2.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.3
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.1.4
mawk -v x=0.7 '$5 > x' B_1.3.3.2.3.2.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.1
mawk -v x=0.65 '$5 > x' B_1.3.3.2.3.2.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.2
mawk -v x=0.6 '$5 > x' B_1.3.3.2.3.2.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.3
mawk -v x=0.5 '$5 > x' B_1.3.3.2.3.2.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.2.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.2.2.4
mawk -v x=0.7 '$5 > x' B_1.3.3.2.3.4.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.1.1
mawk -v x=0.65 '$5 > x' B_1.3.3.2.3.4.1.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.1.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.1.2
mawk -v x=0.7 '$5 > x' B_1.3.3.2.3.4.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.2.1
mawk -v x=0.65 '$5 > x' B.1.3.3.2.3.4.2.imiss | cut -f1 > lowDP.indv
vcftools --vcf B.1.3.3.2.3.4.2.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out B.1.3.3.2.3.4.2.2

SOL.filter2.sh
