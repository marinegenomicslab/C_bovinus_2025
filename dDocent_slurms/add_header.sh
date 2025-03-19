#!/bin/bash

i=$1

echo "Checking $i"

if [ -a $i.bam ] ; then
    if [ $(samtools view -H $i.bam | head -n1 | grep HD | wc -l) -eq 0 ]; then
    echo "Adding HD heading to $i"
    samtools view -h $i.bam | sed -e '1 s/^/@HD\tVN:1.4\tSO:coordinate\n/' | samtools view -b > ${i}_tmp.bam; mv ${i}_tmp.bam ${i}.bam; fi

elif [ -a ${i}-RG.bam ] ; then
    if [ $(samtools view -H $i-RG.bam | head -n1 | grep HD | wc -l) -eq 0 ]; then 
    echo "Adding HD heading to $i"
    samtools view -h $i-RG.bam | sed -e '1 s/^/@HD\tVN:1.4\tSO:coordinate\n/' | samtools view -b > ${i}.bam; fi

else echo "bam file for $i not found"; fi
