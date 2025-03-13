###### Processing the Data from the Cyprinodon Library 1 preparation ######
### Workstations ###
#Move to the directory
```{bash}```
cd ~/Workspace/Robert/Cyprinodon

#Make directories
```{bash}```
mkdir data mapping reference

#Demultiplexing
```{bash}```
cd data
mkdir I2 I4 I7 I10
nano demultiplex_Cyp_Lib1.txt
# Paste in sample name \t index sequence \t barcode sequence
demultiplex demultiplex.pl -i demultiplex_Cyp_Lib1.txt -o Extract_Lib1 -p . -d /home/DATA/GAMBUSIA/Cyp_Lib1
# hand removed the 2nd enzyme
# hand copied each index into an extraction file in each directory

#Demultiplex samples
```{bash}```
ls -d I* | while read i; do
cd $i
sh Extract_Lib1 &
cd .. & done

#Make a record of the samples extracted and how many reads in each
```{bash}```
if [ -a read_log.txt ]; then rm read_log.txt; fi
ls -d I* | while read i; do
DIR=$(echo $i | cut -f1 -d"/")
IND=$(ls $i/*.log | cut -f1 -d"_")
BAR=$(ls $i/*_barcodes.txt | cut -f1 -d"_" | cut -f2 -d"/")
awk -v DIR=$DIR -v IND=$BAR -v FS=" " -v OFS="\t" 'NF==5{print IND, $0, DIR}' $i/*.log | tail -n+3 >> read_log.txt
done
sed -i '1 s/^/Index\tBarcode\tTotal\tNo_Radtag\tLow_Quality\tRetained\tDirectory\n/' read_log.txt

#Adding sample names to read_log.txt
```{R}```
log.dat <- as.data.frame(read.table("read_log.txt", head=T, sep="\t"))
meta.dat <- read.table("demultiplex_Cyp_Lib1.txt", head=F)
tmp.v <- NULL
for(i in 1:nrow(log.dat)){
tmp.name <- meta.dat[which(meta.dat[,2] == log.dat$Barcode & meta.dat[,3] == log.dat$Index), 1]
tmp.v <- c(tmp.v, tmp.name)
}
log.dat$Sample <- tmp.v

#Trimming reads (Set to 10 threads for each)
```{bash}```
screen -S Cyprinodon
ls -d I* | while read i; do
cd $i
cp ~/bin/trim_config.file .
dDocent trim_config.file & done

#Making subdirectories in the reference directory
```{bash}```
cd ../reference
mkdir K1_1_K2_2 K1_2_K2_1 K1_2_K2_2 ref_set test_set

#Add a list of files to put in the reference set
#Look at read_log.txt to find those with the most reads to put into the reference set
#Look at read_log.txt to find good samples for the test set
```{bash}```
nano ref_set.txt
nano test_set.txt

#Copy files to the directories
```{bash}```
ls *_set.txt | while read i; do
echo "Processing $i"
SET=$(echo $i | cut -f1 -d".")
cd $SET
cat ../$i | while read j; do
cp ../../data/I*/${j}*.fq.gz .
done
cd ..
done

#Making the reference directories
```{bash}```
ls -d K1* | while read i; do
echo "making" $i
K1=$(echo $i | sed -e 's/_/\t/g' | cut -f2)
K2=$(echo $i | sed -e 's/_/\t/g' | cut -f4)
cd $i
for j in $(seq 0.8 0.02 0.98); do
mkdir c_${j}
cd c_${j}
cp -s ../../ref_set/*.fq.gz .
sed "s/0.80/$j/g" $WORK/bin/slurm/ref_config.file | awk -v var1="$K1" -v var2="$K2" '$0 ~ /K1/{print $0; getline; gsub($0,var1)}$0 ~ /K2/{print $0; getline; gsub($0,var2)}{print $0}' > ref_config.file
cd ..
done
cd ..
done

#Config file check
grep -rn -A1 Clustering_Similarity K1_1_K2_2/c_0.*/ref_config.file

#For making the reference config files if they were not right in the last loop
##Not used unless making the directories goes wrong##
ls -d K1* | while read i; do
K1=$(echo $i | sed -e 's/_/\t/g' | cut -f2)
K2=$(echo $i | sed -e 's/_/\t/g' | cut -f4)
for j in $(seq 0.8 0.02 0.98); do
sed "s/0.80/$j/g" $WORK/bin/slurm/ref_config.file | awk -v var1="$K1" -v var2="$K2" '$0 ~ /K1/{print $0; getline; gsub($0,var1)}$0 ~ /K2/{print $0; getline; gsub($0,var2)}{print $0}' > ${i}/c_${j}/ref_config.file
done
done

#Removing the c-value directories if necessary
##Not used unless making the directories goes wrong##
ls -d K1* | while read i; do 
cd $i
rm -r *
cd ..
done

#Running the reference building
ls -d K1* | while read i; do
cd $i
for j in $(seq 0.8 0.02 0.98); do
cd c_${j}
while [ $(squeue | grep afields | awk '$5=="R" {print $0} $5=="PD" {print $0}' | wc -l) -ge 6 ]; do sleep 60; done;
echo "starting" $i $j
sbatch -p jgoldq,normal $WORK/bin/slurm/dDocent_ref_assembly.slurm
sleep 2
cd ..; done
cd ..; done

#Checking if references built
ls -d K1_*_K2_*/c_* | while read i; do if [ ! -s $i/reference.fasta ]; then echo "$i missing reference"; fi; done

#Running on ones without a reference
ls -d K1_*_K2_*/c_* | while read i; do 
if [ ! -s $i/reference.fasta ]; then
echo "restarting $i"
cd $i
sbatch -p jgold,normal $WORK/bin/slurm/dDocent_ref_assembly.slurm
cd ../..
sleep 2
fi
done
#Done by hand due to node issues
sbatch -w hpcc07 $WORK/bin/slurm/dDocent_ref_assembly.slurm

#Removing extraneous files from reference making and preparing for mapping
mkdir tmp
ls -d K1* | while read i; do
cd $i
for j in $(seq 0.80 0.02 0.98); do
if [ ! -s c_$j/reference.fasta ]; then echo "${i}/${j} missing reference"; continue; fi
echo "starting" $i $j
cd c_${j}
mv reference.fasta* ../../tmp
mv ref_config.file ../../tmp
rm -fr *
mv ../../tmp/* .
awk '$0 ~ /Assembly/{print $0; getline; gsub($0,"no")}$0 ~ /Reads/{print $0; getline; gsub($0,"yes")}{print $0}' ref_config.file > map_config.file
cp -s ../../test_set/*.fq.gz .
cd ..; done
cd ..; done
rm -r tmp

#Running the reference testing
ls -d K1* | while read i; do
for j in $(seq 0.80 0.02 0.98); do
if [ ! -s ${i}/c_${j}/reference.fasta ]; then echo "Reference missing from $i $j"; continue; fi
cd ${i}/c_${j}/
while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done
echo "starting" $i $j
sbatch -p jgold,normal $WORK/bin/slurm/dDocent_mapping_w_stats.slurm
done
cd ../..
done
#Done by hand due to node issues
sbatch -w hpcc07 $WORK/bin/slurm/dDocent_mapping_w_stats.slurm

#Gathering mapping stats
echo -e "K_values\tc-value\tFile\tReference\tTotal\tQC\tMapped\tPaired"> mapping_summary_all.txt
ls -d K1* | while read i; do
for j in $(seq 0.80 0.02 0.98); do
echo $i $j
if [ ! -s ${i}/c_${j}/reference.fasta ]; then echo "Reference missing from $i $j"; continue; fi
cd ${i}/c_${j}/
tail -n+2 bwa_mapping_summary.txt | awk -v Kval=$i -v Cval=$j '{print Kval, Cval, $0}'  >> ../../mapping_summary_all.txt
cd ../..;
done
done

tail -n+2 mapping_summary_all.txt |cut -f1 -d" " | awk 'BEGIN{FS="_"; OFS="\t"}{print $2, $4}' > K_values.txt
tail -n+2 mapping_summary_all.txt | cut -d" " -f2- | sed 's/ \+ / /g' | sed 's/ /\t/g' > data.txt
paste -d"\t" K_values.txt data.txt > mapping_summary_df.txt
sed -i  "1 s/^/K1\tK2\tcvalue\tFile\tTotal\tQC\tMapped\tPaired\n/" mapping_summary_df.txt

#Gathering reference stats
echo -e "K-values\tcvalue\tcontigs"> ref_summary.txt
ls -d K1* | while read i; do
for j in $(seq 0.80 0.02 0.98); do
echo $i $j
if [ ! -s ${i}/c_${j}/reference.fasta ]; then echo "Reference missing from $i $j"; continue; fi
cd ${i}/c_${j}/
TMP=$(grep ">" reference.fasta | wc -l)
echo -e "$i\t$j\t$TMP" >> ../../ref_summary.txt
cd ../..
done
done

echo -e "K1\tK2\tc-value\tcontigs"> ref_summary_df.txt
tail -n+2 ref_summary.txt | awk 'BEGIN{FS="_|\t"; OFS="\t"}{print $2, $4, $5, $6}' >> ref_summary_df.txt

#Transferring data to Earth for analysis
scp ref_summary_df.txt.txt mapping_summary_df.txt ref_summary.txt ref_summary_df.txt afields@earth.ad.tamucc.edu:/home/afields/Workspace/Robert/Cyprinodon/reference

```{R}```
#define function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} #http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

#Import libraries
library('tidyr')
library('dplyr')
library('ggplot2')

#Import and combine data
dat.map <- read.table("mapping_summary_df.txt", head=T)
dat.con <- read.table("ref_summary_df.txt", head=T)
names(dat.con)[3] <- "cvalue"
dat <- merge(dat.map, dat.con, by=c("K1","K2", "cvalue"))
head(dat.map)
head(dat.con)
head(dat)

#Manipulate data
dat$QC_per <- dat$QC/(dat$Total*2)
dat$Mapped_per <- dat$Mapped/(dat$Total*2)
dat$Paired_per <- dat$Paired/(dat$Total*2)
head(dat)

#Plot number of reads for each value of c for each K1 and K2 combination
ggplot(dat, aes(x = cvalue, y = contigs)) +
  theme(panel.grid.major=element_line(color="grey"), 
        panel.border=element_rect(color="black", fill=NA), 
        axis.text.x=element_text(angle=90,hjust=1)) +
  geom_point() +
  geom_line() +
  facet_grid(K1 ~ K2, labeller=label_both) +
  labs(x = "% similarity c", y = "no. of contigs in reference")

#histogram of the number of contigs per reference
hist(dat.con$contigs, breaks=10, xlab="Number of contigs in references", main="Histogram of the number of references in each combination")
abline(v=mean(dat.con$contigs), col="red4")
abline(v=median(dat.con$contigs), col="mediumblue")
mtext("mean", 3, adj=0.95, col="red4")
mtext("median", 3, line=-1, adj=0.95, col="mediumblue")

#Boxplot of number of contigs per contig vs cvalue
boxplot(dat$contigs ~ dat$cvalue)
box<-boxplot(dat$contigs ~ dat$cvalue, plot=F)
box

#Looking at the quartiles
quart<-cbind(as.numeric(noquote(box$names)), box$stat[1,], box$stat[3,], box$stat[5,])
colnames(quart)<-c("c","low", "med", "high")
quart
plot(quart[,1],quart[,2], col="blue", pch=16, ylim=c(min(quart[,2:4]),max(quart[,2:4])))
points(quart[,1],quart[,3], col="orange", pch=16)
points(quart[,1],quart[,4], col="red", pch=16)

#Computing the lag based upon 0.2 difference
lag1<-rbind(diff(quart[,4], lag=1), diff(quart[,3], lag=1), diff(quart[,2], lag=1))
rownames(lag1)<-c("75%", "Median", "25%")
colnames(lag1)<-c("0.80-0.82", "0.82-0.84", "0.84-0.86", "0.86-0.88", "0.88-0.90", "0.90-0.92", "0.92-0.94", "0.94-0.96", "0.96-0.98")
lag1

#Computing the lag based upon 0.4 difference
lag2<-rbind(diff(lag1[1,], lag1=1), diff(lag1[2,], lag1=1), diff(lag1[3,], lag1=1))
rownames(lag2)<-rownames(lag1)
colnames(lag2)<-c("0.80-0.84", "0.82-0.86", "0.84-0.88", "0.86-0.90", "0.88-0.92", "0.90-0.94", "0.92-0.96", "0.94-0.98")
lag2

par(mfrow = c(2, 1))

plot(lag1[1,], col="red", pch=16, ylim=c(min(lag1),max(lag1)), xaxt="n", xlab=NULL)
axis(1, at= 1:9, labels= colnames(lag1), las=2)
points(lag1[2,], col="orange", pch=16)
points(lag1[3,], col="blue", pch=16)
legend(x=1, y=max(lag1[1,]), legend=c("75% Q", "50% Q", "25% Q"), pch=16, col=c("red", "orange", "blue"))

plot(lag2[1,], col="red", pch=16, ylim=c(min(lag2),max(lag2)), xaxt="n", xlab=NULL)
axis(1, at= 1:8, labels= colnames(lag2), las=2)
points(lag2[2,], col="orange", pch=16)
points(lag2[3,], col="blue", pch=16)
legend(x=1, y=max(lag2[1,]), legend=c("75% Q", "50% Q", "25% Q"), pch=16, col=c("red", "orange", "blue"))

tiff("Lag_graphics.tif", res=200, height=4000, width=2000)
par(mfrow = c(2, 1))

plot(lag1[1,], col="red", pch=16, ylim=c(min(lag1),max(lag1)), xaxt="n", xlab=NULL)
axis(1, at= 1:9, labels= colnames(lag1), las=2)
points(lag1[2,], col="orange", pch=16)
points(lag1[3,], col="blue", pch=16)
legend(x=1, y=max(lag1[1,]), legend=c("75% Q", "50% Q", "25% Q"), pch=16, col=c("red", "orange", "blue"))

plot(lag2[1,], col="red", pch=16, ylim=c(min(lag2),max(lag2)), xaxt="n", xlab=NULL)
axis(1, at= 1:8, labels= colnames(lag2), las=2)
points(lag2[2,], col="orange", pch=16)
points(lag2[3,], col="blue", pch=16)
legend(x=1, y=max(lag2[1,]), legend=c("75% Q", "50% Q", "25% Q"), pch=16, col=c("red", "orange", "blue"))
dev.off()

#Mapping Comparisons
dat.map$Map_Per <- dat.map$Mapped/(dat.map$Total*2)
dat.map$Pair_Per <- dat.map$Paired/(dat.map$Total*2)

par(mfrow = c(3, 1))
boxplot(Map_Per ~ cvalue, data=dat.map[which(dat.map$K1 == 1 & dat.map$K2 == 2),], ylim=c(0.8,1))
boxplot(Map_Per ~ cvalue, data=dat.map[which(dat.map$K1 == 2 & dat.map$K2 == 1),], ylim=c(0.8,1))
boxplot(Map_Per ~ cvalue, data=dat.map[which(dat.map$K1 == 2 & dat.map$K2 == 2),], ylim=c(0.8,1))

par(mfrow = c(3, 1))
hist(dat.map$Map_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 1 & dat.map$K2 == 2)], col="red4", breaks=seq(0,1,0.01), xlab="Percent Mapped", main="K1 = 1 and K2 = 2", xlim=c(0,1))
hist(dat.map$Map_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 1)], col="darkgreen", breaks=seq(0,1,0.01), xlab="Percent Mapped", main="K1 = 2 and K2 = 1", xlim=c(0,1))
hist(dat.map$Map_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 2)], col="mediumblue", breaks=seq(0,1,0.01), xlab="Percent Mapped", main="K1 = 2 and K2 = 2", xlim=c(0,1))


which(dat.map$Map_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 1 & dat.map$K2 == 2)] < 0.5)

tiff("paired_dist.tif", res=200, height=6000, width=2000)
par(mfrow = c(3, 1))
hist(dat.map$Pair_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 1 & dat.map$K2 == 2)], col="red4", breaks=seq(0,1,0.01), xlab="Percent Properly Paired", main="K1 = 1 and K2 = 2", xlim=c(0,1))
abline(v=mean(dat.map$Pair_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 1 & dat.map$K2 == 2)]), col="black")
abline(v=median(dat.map$Pair_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 1 & dat.map$K2 == 2)]), col="black", lty=2)
hist(dat.map$Pair_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 1)], col="darkgreen", breaks=seq(0,1,0.01), xlab="Percent Properly Paired", main="K1 = 2 and K2 = 1", xlim=c(0,1))
abline(v=mean(dat.map$Pair_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 1)]), col="black")
abline(v=median(dat.map$Pair_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 1)]), col="black", lty=2)
hist(dat.map$Pair_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 2)], col="mediumblue", breaks=seq(0,1,0.01), xlab="Percent Properly Paired", main="K1 = 2 and K2 = 2", xlim=c(0,1))
abline(v=mean(dat.map$Pair_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 2)]), col="black")
abline(v=median(dat.map$Pair_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 2)]), col="black", lty=2)
dev.off()

#### Reference assembly complete ####
#### HPC ####
cd $WORK/Workspace/Robert/Cyprinodon/reference
cp K1_1_K2_2/c_0.88/reference.fasta ./Cyprinodon_reference.fasta

#Mapping
#Mapping data
#Run on normal nodes, scripted to limit to 8 nodes
cd ../data
find ./I* -type d > ../mapping/dirs.txt; cd ../mapping; xargs mkdir -p < dirs.txt; rm dirs.txt

ls | while read i; do cd $i; cp -s ../../data/$i/*.fq.gz .; cp -s ../../reference/Cyprinodon_reference.fasta ./reference.fasta; cp -s $WORK/bin/slurm/map_config.file .; cd ..; done
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/bin/slurm/dDocent_mapping.slurm; cd ..; sleep 5; done

#Filtering bam files
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/bin/slurm/dDocent_bamfiltering.slurm; cd ..; sleep 5; done

#Making coverage files for all data
#Run on normal nodes, scripted to limit to 10 nodes
ls | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/bin/slurm/dDocent_cov.slurm; cd ..; sleep 5; done

#Preparing folder for SNP calling
mkdir ../SNP_calling
cd ../SNP_calling
mkdir tmp
ls ../mapping | while read i; do echo $i; cp -s ../mapping/$i/*.bam* .; rm *cat*; done
ls ../mapping | while read i; do echo $i; cp -s ../mapping/$i/*.cov.stats .; done
ls ../mapping | while read i; do echo $i; cat ../mapping/$i/mapped.bed >> all_mapped.bed; done
sbatch $WORK/bin/slurm/dDocent_bed.slurm

ls *.bam | sed 's/.bam//g' > namelist
cut -f1 -d "_" namelist > p; paste namelist p > popmap; rm p
sbatch $WORK/bin/slurm/dDocent_cat.slurm

#Splitting cat.bam for SNP calling
#Number of nodes used is determined by this split script
sbatch -p jgoldq,normal --export=NODES=8 $WORK/bin/slurm/dDocent_split.slurm
for i in $(seq 1 8); do cd $i.node; cp -fs ../cat-RRG.bam* .; cd ..; done
ls -d *.node | while read i; do cd $i; cp -s ../../reference/Cyprinodon_reference.fasta ./reference.fasta; cd .. ; done

#SNP calling on normal nodes, scripted to limit to 8 nodes
ulimit -s 81920
ls -d *.node | while read i; do while [ $(squeue | grep afields | awk '$5=="R" {print $0} $5=="PD" {print $0}' | wc -l) -ge 8 ]; do sleep 60; done; cd $i; ulimit -s 81920; sbatch -p jgoldq,normal $WORK/bin/slurm/dDocent_freebayes.slurm; cd ..; sleep 2; done

#Checking to see if nodes have the correct number of vcf files with data in them
ls -d *.node | while read i; do echo "checking $i"; cd $i; BEDS=$(ls mapped.*.bed | wc -l);	VCF=$(find . -name "*.vcf" -size +62k | wc -l); if [ $VCF -lt $BEDS ]; then echo $i "did not complete all the vcf files properly"; fi; if [ $(find . -name "*.vcf" | wc -l) -gt $BEDS ]; then echo $i "has too many vcf files present"; fi; cd ..; done

#Checking to see if vcf file have the same number of contigs in the header as the reference.fasta has
ls -d *.node | while read i; do echo "checking $i"; cd $i; ls raw.*.vcf | while read j; do VCF=$(head -n1 $j| grep "##" | wc -l); if [ $VCF -eq 0 ]; then echo $i $j "missing complete header"; echo ${i}/${j} >> ../bad.vcf; fi; done; cd ..; done

#Combine all vcfs in each node
ls -d *.node | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch -p jgoldq,normal $WORK/bin/slurm/dDocent_combine_node.slurm; cd ..; sleep 5; done

while [ $(squeue -j 160579 | wc -l) -gt 1 ]; do sleep 30; done; ls -d *.node | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch -p jgoldq,normal $WORK/bin/slurm/dDocent_combine_node.slurm; cd ..; sleep 5; done

#Preparing to combine all the vcf files
mkdir vcf
cd vcf

ls -d ../*.node | while read i; do 
NODE=$(echo $i | sed 's:../::g' | cut -f1 -d .)
cp -fs $i/cat.vcf ./raw.$NODE.vcf
done

#Combining all the vcf files
sbatch -p jgoldq,normal $WORK/bin/slurm/dDocent_combine.slurm

#Transfer
scp ../popmap TotalRawSNPs.vcf afields@earth.ad.tamucc.edu:/home/afields/Workspace/Robert/Cyprinodon/filtering

#When adding the two samples from the Gambusia library
scp ../popmap TotalRawSNPs.vcf afields@earth.ad.tamucc.edu:/home/afields/Workspace/Robert/Cyprinodon/filtering_v2





















#Copy directories from Earth
scp -r afields@earth.ad.tamucc.edu:/home/afields/Workspace/Robert/Cyprinodon .
rm data/I*/*.rem.*

#Trimming data
cd data
ls -d I* | while read i; do cd $i; cp -s $WORK/bin/slurm/trim_config.file .; cd ..; done
ls -d I* | while read i; do while [ $(squeue | grep afields | wc -l) -ge 8 ]; do sleep 60; done; cd $i; sbatch $WORK/bin/slurm/dDocent_trimming.slurm; cd ..; sleep 2; done

#Making subdirectories in the reference directory
```{bash}```
cd ../reference
mkdir K1_1_K2_2 K1_2_K2_1 K1_2_K2_2 ref_set test_set

#Making the reference directories
ls -d K1* | while read i; do
echo "making" $i
K1=$(echo $i | sed -e 's/_/\t/g' | cut -f2)
K2=$(echo $i | sed -e 's/_/\t/g' | cut -f4)
cd $i
for j in $(seq 0.8 0.02 0.98); do
mkdir c_${j}
cd c_${j}
cp -s ../../ref_set/*.fq.gz .
sed "s/0.80/$j/g" $WORK/bin/slurm/ref_config.file | awk -v var1="$K1" -v var2="$K2" '$0 ~ /K1/{print $0; getline; gsub($0,var1)}$0 ~ /K2/{print $0; getline; gsub($0,var2)}{print $0}' > ref_config.file
cd ..
done
cd ..
done

#For making the reference config files if they were not right in the last loop
##Not used unless making the directories goes wrong##
ls -d K1* | while read i; do
K1=$(echo $i | sed -e 's/_/\t/g' | cut -f2)
K2=$(echo $i | sed -e 's/_/\t/g' | cut -f4)
for j in $(seq 0.8 0.02 0.98); do
sed "s/0.80/$j/g" $WORK/bin/slurm/ref_config.file | awk -v var1="$K1" -v var2="$K2" '$0 ~ /K1/{print $0; getline; gsub($0,var1)}$0 ~ /K2/{print $0; getline; gsub($0,var2)}{print $0}' > ${i}/c_${j}/ref_config.file
done
done

#Removing the c-value directories if necessary
##Not used unless making the directories goes wrong##
ls -d K1* | while read i; do 
cd $i
rm -r *
cd ..
done

#Running the reference building
ls -d K1* | while read i; do
cd $i
for j in $(seq 0.8 0.02 0.98); do
cd c_${j}
while [ $(squeue | grep afields | awk '$5=="R" {print $0} $5=="PD" {print $0}' | wc -l) -ge 8 ]; do sleep 60; done;
echo "starting" $i $j
if [ $(squeue | grep jgoldq | wc -l) -eq 0 ]; then sbatch -p jgoldq $WORK/bin/slurm/dDocent_ref_assembly.slurm;
elif [ $(squeue | grep jgoldq | wc -l) -eq 1 ]; then sbatch $WORK/bin/slurm/dDocent_ref_assembly.slurm; fi
sleep 2
cd ..; done
cd ..; done

#Checking if references built
ls -d K1_*_K2_*/c_* | while read i; do if [ ! -s $i/reference.fasta ]; then echo "$i missing reference"; fi; done

#Running on ones without a reference
ls -d K1_*_K2_*/c_* | while read i; do 
if [ ! -s $i/reference.fasta ]; then
while [ $(squeue | grep afields | awk '$5=="R" {print $0} $5=="PD" {print $0}' | wc -l) -ge 8 ]; do sleep 60; done;
echo "restarting $i"
cd $i
if [ $(squeue | grep jgoldq | wc -l) -eq 0 ]; then sbatch -p jgoldq $WORK/bin/slurm/dDocent_ref_assembly.slurm
elif [ $(squeue | grep jgoldq | wc -l) -eq 1 ]; then sbatch $WORK/bin/slurm/dDocent_ref_assembly.slurm
fi
cd -
sleep 2
fi
done

#Removing extraneous files from reference making and preparing for mapping
mkdir tmp
ls -d K1* | while read i; do
echo "Processing $i"
cd $i
for j in $(seq 0.80 0.02 0.98); do
if [ ! -s c_$j/reference.fasta ]; then continue; fi
cd c_${j}
mv reference.fasta* ../../tmp
mv ref_config.file ../../tmp
rm -fr *
mv ../../tmp/* .
awk '$0 ~ /Assembly/{print $0; getline; gsub($0,"no")}$0 ~ /Reads/{print $0; getline; gsub($0,"yes")}{print $0}' ref_config.file > map_config.file
cp -s ../../test_set/*.fq.gz .
cd ..; done
cd ..; done
rm -r tmp

#Running the reference testing
ls -d K1* | while read i; do
for j in $(seq 0.80 0.02 0.98); do
if [ ! -s ${i}/c_${j}/reference.fasta ]; then echo "Reference missing from $i $j"; continue; fi
cd ${i}/c_${j}/
while [ $(squeue | grep afields | awk '$5=="R" {print $0} $5=="PD" {print $0}' | wc -l) -ge 8 ]; do sleep 60; done;
echo "starting" $i $j
sbatch -p jgoldq,normal $WORK/bin/slurm/dDocent_mapping_w_stats.slurm
sleep 2
cd ../..;
done
done

#Gathering mapping stats
echo -e "K_values\tc-value\tFile\tReference\tTotal\tQC\tMapped\tPaired"> mapping_summary_all.txt
ls -d K1* | while read i; do
for j in $(seq 0.80 0.02 0.98); do
echo $i $j
if [ ! -s ${i}/c_${j}/reference.fasta ]; then echo "Reference missing from $i $j"; continue; fi
cd ${i}/c_${j}/
tail -n+2 bwa_mapping_summary.txt | awk -v Kval=$i -v Cval=$j '{print Kval, Cval, $0}'  >> ../../mapping_summary_all.txt
cd ../..;
done
done

tail -n+2 mapping_summary_all.txt |cut -f1 -d" " | awk 'BEGIN{FS="_"; OFS="\t"}{print $2, $4}' > K_values.txt
tail -n+2 mapping_summary_all.txt | cut -d" " -f2- | sed 's/ \+ / /g' | sed 's/ /\t/g' > data.txt
paste -d"\t" K_values.txt data.txt > mapping_summary_df.txt
sed -i  "1 s/^/K1\tK2\tcvalue\tFile\tTotal\tQC\tMapped\tPaired\n/" mapping_summary_df.txt

#Gathering reference stats
echo -e "K-values\tcvalue\tcontigs"> ref_summary.txt
ls -d K1* | while read i; do
for j in $(seq 0.80 0.02 0.98); do
echo $i $j
if [ ! -s ${i}/c_${j}/reference.fasta ]; then echo "Reference missing from $i $j"; continue; fi
cd ${i}/c_${j}/
TMP=$(grep ">" reference.fasta | wc -l)
echo -e "$i\t$j\t$TMP" >> ../../ref_summary.txt
cd ../..
done
done

echo -e "K1\tK2\tc-value\tcontigs"> ref_summary_df.txt
tail -n+2 ref_summary.txt | awk 'BEGIN{FS="_|\t"; OFS="\t"}{print $2, $4, $5, $6}' >> ref_summary_df.txt

scp *_df.txt afields@earth.ad.tamucc.edu:/home/afields/Workspace/Robert/Cyprinodon
