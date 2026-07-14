# Perform mapping to reference genome using BWA
## *Primula grandis* samples
We run them as an array.
```sh
#!/bin/bash
#SBATCH --job-name=EMC_bwa
#SBATCH -o ./logs/bwa.%j.out
#SBATCH -e ./logs/bwa.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=3800
#SBATCH --nice=1000
#SBATCH --partition=main
#SBATCH --array=1-10%10

INDEX=/home/ubuntu/out_JBAT_review4.FINAL_rename.fasta
FORWARD=$(ls /home/ubuntu/trimmed_reads/p_grandis/grandis_*_R1_paired.fastq.gz | head -n $SLURM_ARRAY_TASK_ID | tail -1)
REVERSE=$(echo $FORWARD | sed 's/_R1/_R2/g')
NAME=$(basename $FORWARD _R1_paired.fastq.gz)
OUTPUT=./${NAME}.bam
SAMTOOLS=/home/ubuntu/executables/samtools-1.10/samtools
BWA=/home/ubuntu/executables/bwa-0.7.17/bwa

echo "Starting aligment $NAME with bwa"

$BWA mem -t $SLURM_CPUS_PER_TAKS $INDEX $FORWARD $REVERSE | $SAMTOOLS view -@ $SLURM_CPUS_PER_TASK -b | $SAMTOOLS sort -@ $SLURM_CPUS_PER_TASK -T ${NAME} > $OUTPUT

echo "Finish alignment of $NAME"

$SAMTOOLS index $OUTPUT

echo "Finish indexing of $NAME"
```
## Other *Primula* species
For all the other *Primula* species, we mapped first to subgenome A and then subgenome B. The subgenomes were identified based on SubPhaser.
Run by using following command script as follows to run as a for loop:
```sh
for x in {1..9}; do sbatch bwa.sh $x; done
```
```sh
#!/bin/bash
#SBATCH --job-name=EMC_bwa
#SBATCH -o ./logs/bwa.%j.out
#SBATCH -e ./logs/bwa.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=3800
#SBATCH --nice=1000
#SBATCH --partition=main

SAMTOOLS=/home/ubuntu/executables/samtools-1.10/samtools
BWA=/home/ubuntu/executables/bwa-0.7.17/bwa

INDEX_SUBGEN_A=/home/ubuntu/p_grandis_ref_genome/subgenomes_fasta/out_JBAT_review4.FINAL_rename_linear_Subgenome_A.fasta
INDEX_SUBGEN_B=/home/ubuntu/p_grandis_ref_genome/subgenomes_fasta/out_JBAT_review4.FINAL_rename_linear_Subgenome_B.fasta
dir_trimmed=/home/ubuntu/trimmed_reads/primula

NUM=${1?Error: no chrom number given}

SAMPLES=$(ls $dir_trimmed/Plate$NUM | grep -oE 'P00'$NUM'_W\w\w\w' | sort -u)

for i in ${SAMPLES[@]}; do
echo "Starting aligment $NAME with bwa"
$BWA mem -t $SLURM_CPUS_PER_TASK $INDEX_SUBGEN_A $dir_trimmed/Plate$NUM/''$i'_R1_paired.fastq.gz' $dir_trimmed/Plate$NUM/''$i'_R2_paired.fastq.gz' | \
$SAMTOOLS view -@ $SLURM_CPUS_PER_TASK -b | $SAMTOOLS sort -@ $SLURM_CPUS_PER_TASK -T $i > $i.SubGenA.bam
echo "Finish alignment of $NAME"
$SAMTOOLS index $i.SubGenA.bam
echo "Finish indexing of $NAME"; done

for i in ${SAMPLES[@]}; do
echo "Starting aligment $NAME with bwa"
$BWA mem -t $SLURM_CPUS_PER_TASK $INDEX_SUBGEN_A $dir_trimmed/Plate$NUM/''$i'_R1_paired.fastq.gz' $dir_trimmed/Plate$NUM/''$i'_R2_paired.fastq.gz' | \
$SAMTOOLS view -@ $SLURM_CPUS_PER_TASK -b | $SAMTOOLS sort -@ $SLURM_CPUS_PER_TASK -T $i > $i.SubGenB.bam
echo "Finish alignment of $NAME"
$SAMTOOLS index $i.SubGenB.bam
echo "Finish indexing of $NAME"; done
```
