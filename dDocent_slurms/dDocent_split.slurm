#!/bin/bash

#SBATCH -J dDocent_split		             # Name of the job
#SBATCH -o dDocent_split.out                 # Name of file that will have program output
#SBATCH -e dDocent_split.err                 # Name of the file that will have job errors, if any
#SBATCH -N 1                                 # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 64                                 # Number of cores ( my test allocated 2 per node )
#SBATCH -p jgoldq,normal                               # Partition
#SBATCH --mail-user=andrew.fields@tamucc.edu
#SBATCH --mail-type=begin                    # email me when the job starts
#SBATCH --mail-type=end                      # email me when the job ends
#SBATCH --time=24:00:00

#module load ddocent/2.9.1

mkdir tmp
cat *.cov.stats | sort -k1,1 -k2,2n --temporary-directory=./tmp | bedtools merge -i - -c 4 -o sum > cov.stats

NUMProc=$(echo "64*$NODES" | bc)
DP=$(mawk '{print $4}' cov.stats | sort -rn | perl -e '$d=.00005;@l=<>;print $l[int($d*@l)]')
CC=$(mawk -v x=$DP '$4 < x' cov.stats | mawk '{len=$3-$2;lc=len*$4;tl=tl+lc} END {OFMT = "%.0f";print tl/"'$NUMProc'"}')

mawk -v x=$DP '$4 < x' cov.stats |sort -V -k1,1 -k2,2 --temporary-directory=./tmp | mawk -v cutoff=$CC 'BEGIN{i=1}
        {len=$3-$2;lc=len*$4;cov = cov + lc
			if ( cov < cutoff) {x="mapped."i".bed";print $1"\t"$2"\t"$3 > x}
			else {i=i+1; x="mapped."i".bed"; print $1"\t"$2"\t"$3 > x; cov=0}
        }'

for i in $(seq 1 $NODES); do mkdir $i.node; START=$(echo "($i-1)*64+1" | bc); END=$(echo "$i*64" | bc); for j in $(seq $START $END); do	mv mapped.$j.bed $i.node/; done; cd $i.node/; cp -s ../cat-RRG.bam* .; cp -s ../reference.fasta* .; cp -s ../popmap .; cd ..; done

rm -r tmp
