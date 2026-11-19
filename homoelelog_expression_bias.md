# Test for subgenome dominance

To test for subgenome dominance, we checked for homoelelog expression bias (HEB) toward one of the subgenomes of P. grandis as follows. First, we aligned transcriptomic data from leaves, flowers, and floral buds to the reference genome using Bowtie2 v2.3.4.3 (Langmead and Salzberg 2012) with the very-sensitive-local option (-D 20 -R 3 -N 0 -L 20 -i S, 1, 0.50). Secondly, the number of aligned reads mapping to each gene was counted using htseq-count v2.0.5 (Putri et al. 2022) to estimate reads per kilobase per million (RPKM). Thirdly, given that we cannot use the likelihood ratio test proposed by Smith et al. (2019) due to the lack of biological and tissue replicates, we established −3 ≤ log2 HEB ≥ 3 (i.e. 8-fold HEB bias toward Subgenome A or B, respectively) as the threshold of HEB significance.

## RNA alignment to reference genome
```sh
INDEX=out_JBAT_review4.FINAL_rename
RNA_DATA=/home/ubuntu/my_data
FORWARD=$(ls $RNA_DATA/*.qc.R1.paired.fastq | head -n $SLURM_ARRAY_TASK_ID | tail -1)
REVERSE=$(echo $FORWARD | sed 's/R1/R2/g')
NAME=$(basename $FORWARD .qc.R1.paired.fastq)
OUTPUT=./${NAME}.bam
BOWTIE=/home/ubuntu/executables/bowtie2-2.3.4.3-linux-x86_64/bowtie2
BWA=/home/ubuntu/executables/programs/bwa-0.7.17/bwa
SAMTOOLS=/home/ubuntu/executables/samtools-1.10/samtools

bowtie2-build /home/ubuntu/p_grandis_ref_genome/out_JBAT_review4.FINAL_rename.fasta out_JBAT_review4.FINAL_rename

echo "Starting aligment $NAME with bwa"

$BOWTIE -p 16 -q --very-sensitive-local -x $INDEX -1 $FORWARD -2 $REVERSE -S ${NAME}.sam

$SAMTOOLS view -@ $SLURM_CPUS_PER_TASK -b ${NAME}.sam | $SAMTOOLS sort -@ $SLURM_CPUS_PER_TASK -T ${NAME}.sam > $OUTPUT

echo "Finish alignment of $NAME"

$SAMTOOLS index $OUTPUT

echo "Finish indexing of $NAME"

rm ${NAME}.sam
```
## Counts of mapped with htseq-count
```sh
DIR=/home/ubuntu/dominance_subgenomes/smith_2019
GTF=/home/ubuntu/p_grandis_ref_genome/braker_rename.gtf

EXE=/home/ubuntu/executables/progs/miniconda3/envs/emiliano3.7/bin

BAM=(B-grandis_F2 B-grandis_L2 B-grandis_buds1)

for i in ${BAM[@]}; do
$EXE/htseq-count -f bam --order=pos --stranded=no $DIR/bowtie_alignments/20200727.$i.bam $GTF -n $SLURM_CPUS_PER_TASK --max-reads-in-buffer 9223372036854775807 \
--mode=intersection-nonempty -t CDS -i transcript_id -a 40 > $i.txt; done
```
## Check gene lenghts
```sh
BAM=(B-grandis_F2 B-grandis_L2 B-grandis_buds1)
rm *txt
for i in ${BAM[@]}; do
GENES=$(awk '{print $1}' $DIR/count_reads/$i.txt)
for j in ${GENES[@]}; do
grep '\b'$j'\b' $GTF | grep CDS | awk '{sum -= $5 - $4} END {print sum}' | awk '{$NF *= -1; print}' >> gene_lenghts.$i.txt;done ;done
```
## Estimation of RPKM
Calculate the RPKM (Reads per kilo base per million mapped reads) for the RNAseq data.
Formula RPKM = numReads / (geneLength/1000 * totalNumReads/1,000,000)
```sh
DIR=/home/ubuntu/emiliano/4th_chp/8_dominance_subgenomes/smith_2019/count_reads

paste $DIR/htseq_outputs/B-grandis_F2.txt $DIR/gene_lenghts/gene_lenghts.B-grandis_F2.txt | awk 'BEGIN{FS=OFS="\t"} {$2=($2==0)?"NA":$2} 1' | head -n -5 > HTseq_gene_lenghts.B-grandis_F2.tx
t
paste $DIR/htseq_outputs/B-grandis_F2.txt $DIR/gene_lenghts/gene_lenghts.B-grandis_F2.txt | head -n -5 > HTseq_gene_lenghts.B-grandis_F2.txt
NORM_READS=$(awk '{sum+=$2/1000000} END{print sum;}' $DIR/htseq_outputs/B-grandis_F2.txt)
NORM_GENE_LENGHT=$(awk '{sum+=$1/1000} END{print sum;}' $DIR/gene_lenghts/gene_lenghts.B-grandis_F2.txt)
awk $2!="NA" '{$4=$3/1000;$5+=$2/10.242;$6=$4*$5;$7=$2/$6} {print}' HTseq_gene_lenghts.B-grandis_F2.txt > RPKM.B-grandis_F2.txt
head -n -5 /home/ubuntu/emiliano/4th_chp/8_dominance_subgenomes/smith_2019/count_reads/htseq_outputs/B-grandis_F2.txt | awk '{sum+=$2;} END{print sum/1000000}' = 10.242

FILES=(B-grandis_F2)
for i in ${FILES[@]}; do
paste $DIR/htseq_outputs/$i.txt $DIR/gene_lenghts/gene_lenghts.$i.txt | head -n -5 > $i.temp.txt
awk '$2!="0" {$4=$3/1000;$5+=$2/10.242;$6=$4*$5;$7=$2/$6} {print}' OFS="\t" $i.temp.txt | \
awk 'BEGIN{FS=OFS="\t"} $2 == 0 { $4=0; $5=0; $6=0; $7=0 } 1' > $i.RPKM.txt
rm $i.temp.txt; done

FILES=(B-grandis_F2)
for i in ${FILES[@]}; do
paste $DIR/htseq_outputs/$i.txt $DIR/gene_lenghts/gene_lenghts.$i.txt | head -n -5 > $i.temp.txt
awk '$2!="0" {$4=$3/1000;$5=1*1.159739;$6=$4*$5;$7=$2/$6} {print}' OFS="\t" $i.temp.txt | \
awk 'BEGIN{FS=OFS="\t"} $2 == 0 { $4=0; $5=0; $6=0; $7=0 } 1' > $i.RPKM.txt
rm $i.temp.txt; done

FILES=(B-grandis_L2)
for i in ${FILES[@]}; do
paste $DIR/htseq_outputs/$i.txt $DIR/gene_lenghts/gene_lenghts.$i.txt | head -n -5 > $i.temp.txt
awk '$2!="0" {$4=$3/1000;$5=1*0.863818;$6=$4*$5;$7=$2/$6} {print}' OFS="\t" $i.temp.txt | \
awk 'BEGIN{FS=OFS="\t"} $2 == 0 { $4=0; $5=0; $6=0; $7=0 } 1' > $i.RPKM.txt
rm $i.temp.txt; done

FILES=(B-grandis_buds1)
for i in ${FILES[@]}; do
paste $DIR/htseq_outputs/$i.txt $DIR/gene_lenghts/gene_lenghts.$i.txt | head -n -5 > $i.temp.txt
awk '$2!="0" {$4=$3/1000;$5=1*0.709235;$6=$4*$5;$7=$2/$6} {print}' OFS="\t" $i.temp.txt | \
awk 'BEGIN{FS=OFS="\t"} $2 == 0 { $4=0; $5=0; $6=0; $7=0 } 1' > $i.RPKM.txt
rm $i.temp.txt; done
```
```sh
DIR_1=/home/ubuntu/emiliano/p_grandis_ref_genome/mcscan_stuff
DIR_2=/home/ubuntu/emiliano/4th_chp/8_dominance_subgenomes/smith_2019/count_reads/RPKM

LINES=$(cat $DIR_1/pgraA.pgraB.anchors | awk '{print $1,$2}' OFS=',')

for i in ${LINES[@]}; do
echo $i > temp.file
GREP_1=$(awk -F ',' '{print $1}' temp.file)
grep '\b'$GREP_1'\b' $DIR_2/B-grandis_F2.RPKM.txt | awk '{print "'$GREP_1'",$1,$7}' OFS='\t' >> B-grandis_F2.A.temp
GREP_2=$(awk -F ',' '{print $2}' temp.file)
grep '\b'$GREP_2'\b' $DIR_2/B-grandis_F2.RPKM.txt | awk '{print "'$GREP_2'",$1,$7}' OFS='\t' >> B-grandis_F2.B.temp
paste B-grandis_F2.A.temp B-grandis_F2.B.temp | awk '{print $2,$3,$5,$6}' OFS='\t' > B-grandis_F2.AB.temp
awk 'BEGIN{FS=OFS="\t"} $2!="0" {$5=$4/$2} {print}' B-grandis_F2.AB.temp | awk 'BEGIN{FS=OFS="\t"} $2 == 0 {$5=0}1' > B-grandis_F2.HEB.AB
rm temp.file B-grandis_F2.AB.temp;done
rm B-grandis_F2.A.temp B-grandis_F2.B.temp

FILES=(B-grandis_F2 B-grandis_L2 B-grandis_buds1)
for j in ${FILES[@]}; do
for i in ${LINES[@]}; do
echo $i > temp.file
GREP_1=$(awk -F ',' '{print $1}' temp.file)
grep '\b'$GREP_1'\b' $DIR_2/$j.RPKM.txt | awk '{print "'$GREP_1'",$1,$7}' OFS='\t' >> $j.A.temp
GREP_2=$(awk -F ',' '{print $2}' temp.file)
grep '\b'$GREP_2'\b' $DIR_2/$j.RPKM.txt | awk '{print "'$GREP_2'",$1,$7}' OFS='\t' >> $j.B.temp
paste $j.A.temp $j.B.temp | awk '{print $2,$3,$5,$6}' OFS='\t' > $j.AB.temp
awk 'BEGIN{FS=OFS="\t"} $2!="0" {$5=$4/$2} {print}' $j.AB.temp | awk 'BEGIN{FS=OFS="\t"} $2 == 0 {$5=0}1' > $j.HEB.AB
rm temp.file $j.AB.temp;done
rm $j.A.temp $j.B.temp; done

DIR_1=/home/ubuntu/emiliano/p_grandis_ref_genome/mcscan_stuff
DIR_2=/home/ubuntu/emiliano/4th_chp/8_dominance_subgenomes/smith_2019/count_reads/htseq_outputs_gene_lenght

LINES=$(cat $DIR_1/pgraA.pgraB.anchors | awk '{print $1,$2}' OFS=',')

for i in ${LINES[@]}; do
echo $i > temp.file
GREP_1=$(awk -F ',' '{print $1}' temp.file)
grep '\b'$GREP_1'\b' $DIR_2/B-grandis_L2.txt | awk '{print "'$GREP_1'",$1,$2,$3}' OFS='\t' >> B-grandis_L2.A.temp
GREP_2=$(awk -F ',' '{print $2}' temp.file)
grep '\b'$GREP_2'\b' $DIR_2/B-grandis_L2.txt | awk '{print "'$GREP_2'",$1,$2,$3}' OFS='\t' >> B-grandis_L2.B.temp
paste B-grandis_L2.A.temp B-grandis_L2.B.temp | awk '{print $1,$3,$4,$5,$7,$8}' OFS='\t' > B-grandis_L2.AB.temp
rm temp.file;done
rm B-grandis_L2.A.temp B-grandis_L2.B.temp

for i in ${LINES[@]}; do
echo $i > temp.file
GREP_1=$(awk -F ',' '{print $1}' temp.file)
grep '\b'$GREP_1'\b' $DIR_2/B-grandis_buds1.txt | awk '{print "'$GREP_1'",$1,$2,$3}' OFS='\t' >> B-grandis_buds1.A.temp
GREP_2=$(awk -F ',' '{print $2}' temp.file)
grep '\b'$GREP_2'\b' $DIR_2/B-grandis_buds1.txt | awk '{print "'$GREP_2'",$1,$2,$3}' OFS='\t' >> B-grandis_buds1.B.temp
paste B-grandis_buds1.A.temp B-grandis_buds1.B.temp | awk '{print $1,$3,$4,$5,$7,$8}' OFS='\t' > B-grandis_buds1.AB.temp
rm temp.file;done
rm B-grandis_buds1.A.temp B-grandis_buds1.B.temp
```
## R plot
```R
library(tidyverse)
library(ggpubr)
library(ggplot2)

setwd("~/r_analysis/primula_grandis/databases/dominance")

#flower <- read.table("B-grandis_F2.MQ30.HEB.AB")
#leaf <- read.table("B-grandis_L2.MQ30.HEB.AB")
#bud <- read.table("B-grandis_buds1.MQ30.HEB.AB")

flower <- read.table("B-grandis_F2.HEB.AB")
flower <- flower %>% mutate(log_2 = log2(V5))
flower <- flower %>% 
  filter_all(all_vars(!is.infinite(.)))
                                                     
leaf <- read.table("B-grandis_L2.HEB.AB")
leaf <- leaf %>% mutate(log_2 = log2(V5))
leaf <- leaf %>% 
  filter_all(all_vars(!is.infinite(.)))

bud <- read.table("B-grandis_buds1.HEB.AB")
bud <- bud %>% mutate(log_2 = log2(V5))
bud <- bud %>% 
  filter_all(all_vars(!is.infinite(.)))

BINWIDE <- 1
hist_flowers <- hist(flower$log_2, breaks = seq(-20, 20, by = BINWIDE),
     col = "gray", xlab = "HEBS: Positive=Shift Towards Other in Leaf", ylab = "Count", xlim = c(-10, 10))
axis(1, at = c(-10,-8,-6,-4,-2,0,2,4,6,8,10))

hist_leaf <- hist(leaf$log_2, breaks = seq(-20, 20, by = BINWIDE),
     col = "gray", xlab = "HEBS: Positive=Shift Towards Other in Leaf", ylab = "Count", xlim = c(-10, 10))
axis(1, at = c(-10,-8,-6,-4,-2,0,2,4,6,8,10))

hist_bud <- hist(bud$log_2, breaks = seq(-20, 20, by = BINWIDE),
     col = "gray", xlab = "HEBS: Positive=Shift Towards Other in Floral Buds", ylab = "Count", xlim = c(-10, 10),)
axis(1, at = c(-10,-8,-6,-4,-2,0,2,4,6,8,10))

sum(flower$log_2 >= 3,na.rm=TRUE)
sum(flower$log_2 <= -3,na.rm=TRUE)

sum(leaf$log_2 >= 3,na.rm=TRUE)
sum(leaf$log_2 <= -3,na.rm=TRUE)

sum(bud$log_2 >= 3,na.rm=TRUE)
sum(bud$log_2 <= -3,na.rm=TRUE)
```
