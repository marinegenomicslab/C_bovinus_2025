###### Making a reference for C. bovinus ######
#Trimming data
```{bash}```
cd data
ls -d I* | while read i; do cd $i; cp -s $WORK/slurm/trim_config.file .; cd ..; done
ls -d I* | while read i; do cd $i; sbatch $WORK/slurm/dDocent_trimming.slurm; cd ..; sleep 2; done

#Making subdirectories in the reference directory
```{bash}```
cd ../reference
mkdir K1_1_K2_2 K1_2_K2_1 K1_2_K2_2 ref_set test_set

#Add a list of files to put in the reference set
#Look at read_log.txt to find those with the most reads to put into the reference set
	# Took the 4 with the most reads from each location as long as they did not have more than 10 million reads
#Look at read_log.txt to find good samples for the test set
	# Took all samples above the average number of reads for G.nobilis
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
		sed "s/0.80/$j/g" $WORK/slurm/ref_config.file | awk -v var1="$K1" -v var2="$K2" '$0 ~ /K1/{print $0; getline; gsub($0,var1)}$0 ~ /K2/{print $0; getline; gsub($0,var2)}{print $0}' > ref_config.file
		cd ..
	done
	cd ..
done

#Config file check
```{bash}```
grep -rn -A1 Clustering_Similarity K1_1_K2_2/c_0.*/ref_config.file

#Running the reference building
```{bash}```
ls -d K1* | while read i; do
	cd $i
	for j in $(seq 0.8 0.02 0.98); do
		cd c_${j}
		echo "starting" $i $j
		sbatch $WORK/slurm/dDocent_ref_assembly.slurm
		sleep 2
		cd ..
  	done
	cd ..
done

#Checking if references built
```{bash}```
ls -d K1_*_K2_*/c_* | while read i; do if [ ! -s $i/reference.fasta ]; then echo "$i missing reference"; fi; done

#Removing extraneous files from reference making and preparing for mapping
```{bash}```
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
		cd ..
  	done
	cd ..
done
rm -r tmp

#Running the reference testing
```{bash}```
ls -d K1* | while read i; do
	for j in $(seq 0.80 0.02 0.98); do
		if [ ! -s ${i}/c_${j}/reference.fasta ]; then echo "Reference missing from $i $j"; continue; fi
		cd ${i}/c_${j}/
		echo "starting" $i $j
		sbatch $WORK/slurm/dDocent_mapping_w_stats.slurm
		cd ../..
		sleep 5
	done
done

#Gathering mapping stats
```{bash}```
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
```{bash}```
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

scp ref_summary.txt mapping_summary_df.txt ref_summary.txt ref_summary_df.txt afields@earth.ad.tamucc.edu:/home/afields/Workspace/Robert/Cbovinus/reference

#Import libraries
```{R}```
library('tidyr')
library('dplyr')
library('ggplot2')

#Import and combine data
```{R}```
dat.map <- read.table("mapping_summary_df.txt", head=T)
dat.con <- read.table("ref_summary_df.txt", head=T)
names(dat.con)[3] <- "cvalue"
dat <- merge(dat.map, dat.con, by=c("K1","K2", "cvalue"))
head(dat.map)
head(dat.con)
head(dat)

#Manipulate data
```{R}```
dat$QC_per <- dat$QC/(dat$Total*2)
dat$Mapped_per <- dat$Mapped/(dat$Total*2)
dat$Paired_per <- dat$Paired/(dat$Total*2)
head(dat)

#Plot number of reads for each value of c for each K1 and K2 combination
```{R}```
ggplot(dat, aes(x = cvalue, y = contigs)) +
  theme(panel.grid.major=element_line(color="grey"), 
        panel.border=element_rect(color="black", fill=NA), 
        axis.text.x=element_text(angle=90,hjust=1)) +
  geom_point() +
  geom_line() +
  facet_grid(K1 ~ K2, labeller=label_both) +
  labs(x = "% similarity c", y = "no. of contigs in reference")

#histogram of the number of contigs per reference
```{R}```
hist(dat.con$contigs, breaks=10, xlab="Number of contigs in references", main="Histogram of the number of references in each combination")
abline(v=mean(dat.con$contigs), col="red4")
abline(v=median(dat.con$contigs), col="mediumblue")
mtext("mean", 3, adj=0.95, col="red4")
mtext("median", 3, line=-1, adj=0.95, col="mediumblue")

#Boxplot of number of contigs per contig vs cvalue
```{R}```
boxplot(dat$contigs ~ dat$cvalue)
box<-boxplot(dat$contigs ~ dat$cvalue, plot=F)
box

#Looking at the quartiles
```{R}```
quart<-cbind(as.numeric(noquote(box$names)), box$stat[1,], box$stat[3,], box$stat[5,])
colnames(quart)<-c("c","low", "med", "high")
quart
plot(quart[,1],quart[,2], col="blue", pch=16, ylim=c(min(quart[,2:4]),max(quart[,2:4])), xlab="c-value", ylab="Number of contigs")
points(quart[,1],quart[,3], col="orange", pch=16)
points(quart[,1],quart[,4], col="red", pch=16)

tiff("Quartile_plot.tif", res=200, width=2000, height=2000)
plot(quart[,1],quart[,2], col="blue", pch=16, ylim=c(min(quart[,2:4]),max(quart[,2:4])), xlab="c-value", ylab="Number of contigs")
points(quart[,1],quart[,3], col="orange", pch=16)
points(quart[,1],quart[,4], col="red", pch=16)
dev.off()

#Computing the lag based upon 0.2 difference
```{R}```
lag1<-rbind(diff(quart[,4], lag=1), diff(quart[,3], lag=1), diff(quart[,2], lag=1))
rownames(lag1)<-c("75%", "Median", "25%")
colnames(lag1)<-c("0.80-0.82", "0.82-0.84", "0.84-0.86", "0.86-0.88", "0.88-0.90", "0.90-0.92", "0.92-0.94", "0.94-0.96", "0.96-0.98")
lag1

#Computing the lag based upon 0.4 difference
```{R}```
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

dat.map$Map_Per <- dat.map$Mapped/(dat.map$Total*2)
dat.map$Pair_Per <- dat.map$Paired/(dat.map$Total*2)

par(mfrow = c(3, 1))
boxplot(Map_Per ~ cvalue, data=dat.map[which(dat.map$K1 == 1 & dat.map$K2 == 2),], ylim=c(0.9,1))
boxplot(Map_Per ~ cvalue, data=dat.map[which(dat.map$K1 == 2 & dat.map$K2 == 1),], ylim=c(0.9,1))
boxplot(Map_Per ~ cvalue, data=dat.map[which(dat.map$K1 == 2 & dat.map$K2 == 2),], ylim=c(0.9,1))

par(mfrow = c(3, 1))
hist(dat.map$Map_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 1 & dat.map$K2 == 2)], col="red4", breaks=seq(0,1,0.01), xlab="Percent Mapped", main="K1 = 1 and K2 = 2", xlim=c(0,1))
hist(dat.map$Map_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 1)], col="darkgreen", breaks=seq(0,1,0.01), xlab="Percent Mapped", main="K1 = 2 and K2 = 1", xlim=c(0,1))
hist(dat.map$Map_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 2)], col="mediumblue", breaks=seq(0,1,0.01), xlab="Percent Mapped", main="K1 = 2 and K2 = 2", xlim=c(0,1))

summary(dat.map$Map_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 1 & dat.map$K2 == 2)])
summary(dat.map$Map_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 1)])
summary(dat.map$Map_Per[which(dat.map$cvalue == 0.88 & dat.map$K1 == 2 & dat.map$K2 == 2)])

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

#K1=1, K2=2, c-value= 0.88
