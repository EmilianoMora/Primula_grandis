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
