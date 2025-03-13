#### Reference assembly complete ####
#### HPC ####
cd references
cp K1_1_K2_2/c_0.88/reference.fasta ./Cbovinus_reference.fasta

#Getting data
cd ../data

#Trimming data
cd data
ls -d I* | while read i; do cd $i; cp -s $WORK/bin/slurm/trim_config.file .; cd ..; done
ls -d I* | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/bin/slurm/dDocent_trimming.slurm; cd ..; sleep 2; done

#Mapping
#Making directories
find ./I* -type d > ../mapping/dirs.txt; cd ../mapping; xargs mkdir -p < dirs.txt; rm dirs.txt

#Mapping reads
ls | while read i; do cd $i; cp -s ../../data/$i/*.fq.gz .; cp -s ../../reference/Cbovinus_reference.fasta ./reference.fasta; cp -s $WORK/bin/slurm/map_config.file .; cd ..; done
ls | while read i; do cd $i; sbatch $WORK/bin/slurm/dDocent_mapping.slurm; cd ..; sleep 2; done

#Filtering bam files
ls | while read i; do cd $i; sbatch $WORK/bin/slurm/dDocent_bamfiltering.slurm; cd ..; sleep 5; done

#Making coverage files for all data
#Run on normal nodes, scripted to limit to 10 nodes
ls | while read i; do cd $i; sbatch $WORK/bin/slurm/dDocent_cov.slurm; cd ..; sleep 5; done

#Preparing folder for SNP calling
mkdir ../SNP_calling
cd ../SNP_calling
mkdir tmp
ls ../mapping | while read i; do echo $i; cp -s ../mapping/$i/*.bam* .; rm cat*; done
ls ../mapping | while read i; do echo $i; cp -s ../mapping/$i/*.cov.stats .; done
ls ../mapping | while read i; do echo $i; cat ../mapping/$i/mapped.bed >> all_mapped.bed; done
sbatch $WORK/bin/slurm/dDocent_bed.slurm

ls *.bam | sed 's/.bam//g' > namelist
cut -f1 -d "_" namelist > p; paste namelist p > popmap; rm p
sbatch $WORK/bin/slurm/dDocent_cat.slurm

#Splitting cat.bam for SNP calling
#Number of nodes used is determined by this split script
sbatch -p jgoldq,normal --export=NODES=8 $WORK/bin/slurm/dDocent_split.slurm
#for i in $(seq 1 8); do cd $i.node; cp -fs ../cat-RRG.bam* .; cd ..; done
ls -d *.node | while read i; do cd $i; cp -s ../../reference/Cbovinus_reference.fasta ./reference.fasta; cd .. ; done

#SNP calling on normal nodes, scripted to limit to 8 nodes
ulimit -s 81920
ls -d *.node | while read i; do cd $i; ulimit -s 81920; sbatch -p jgoldq,normal $WORK/bin/slurm/dDocent_freebayes.slurm; cd ..; sleep 2; done

#Checking to see if nodes have the correct number of vcf files with data in them
ls -d *.node | while read i; do echo "checking $i"; cd $i; BEDS=$(ls mapped.*.bed | wc -l);	VCF=$(find . -name "*.vcf" -size +62k | wc -l); if [ $VCF -lt $BEDS ]; then echo $i "did not complete all the vcf files properly"; fi; if [ $(find . -name "*.vcf" | wc -l) -gt $BEDS ]; then echo $i "has too many vcf files present"; fi; cd ..; done

#Checking to see if vcf file have the same number of contigs in the header as the reference.fasta has
ls -d *.node | while read i; do echo "checking $i"; cd $i; ls raw.*.vcf | while read j; do VCF=$(head -n1 $j| grep "##" | wc -l); if [ $VCF -eq 0 ]; then echo $i $j "missing complete header"; echo ${i}/${j} >> ../bad.vcf; fi; done; cd ..; done

#Combine all vcfs in each node
ls -d *.node | while read i; do cd $i; sbatch $WORK/bin/slurm/dDocent_combine_node.slurm; cd ..; sleep 5; done

#Preparing to combine all the vcf files
mkdir vcf
cd vcf

ls -d ../*.node | while read i; do 
NODE=$(echo $i | sed 's:../::g' | cut -f1 -d .)
cp -fs $i/cat.vcf ./raw.$NODE.vcf
done

#Combining all the vcf files
sbatch $WORK/bin/slurm/dDocent_combine.slurm
