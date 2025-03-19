#!/bin/bash

#SBATCH -J dDocent_cov		             # Name of the job
#SBATCH -o dDocent_cov.out                   # Name of file that will have program output
#SBATCH -e dDocent_cov.err                   # Name of the file that will have job errors, if any
#SBATCH -N 1                                 # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 64                                # Number of cores ( my test allocated 2 per node )
#SBATCH -p normal                               # Partition
#SBATCH --mail-user=andrew.fields@tamucc.edu
#SBATCH --mail-type=begin                    # email me when the job starts
#SBATCH --mail-type=end                      # email me when the job ends
#SBATCH --time=24:00:00

#module load ddocent/2.9.1

ls *-RG.bam > bamlist.list
samtools merge -@20 -b bamlist.list -f cat-RRG.bam
samtools index cat-RRG.bam
wait
bedtools merge -i cat-RRG.bam -bed >  mapped.bed

mawk -v OFS='\t' {'print $1,$2'} reference.fasta.fai > genome.file
cat namelist | parallel -j 10 "bedtools coverage -b {}-RG.bam -a mapped.bed -counts -sorted -g genome.file > {}.cov.stats"
