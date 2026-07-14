# Perform trimming of short read Illumina sequences with Trimmomatic v0.38
## Trimming all *Primula grandis* samples that were located in different sequencing plates
```sh
#!/bin/bash
#SBATCH --job-name=EMC_trimm
#SBATCH -o ./logs/trimm.%j.out
#SBATCH -e ./logs/trimm.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=3800
#SBATCH --nice=1000
#SBATCH --partition=main

TRIMMOMATIC=/home/ubuntu/executables/Trimmomatic-0.38
dir_reads=/home/ubuntu/raw_data/p_grandis
dir_trimmed=/home/ubuntu/trimmed_reads/p_grandis

mkdir unpaired_reads

SAMPLES=(WB10 WC09 WD08 WE08)

for i in ${SAMPLES[@]}; do
java -jar $TRIMMOMATIC/trimmomatic-0.38.jar PE $dir_reads'/P004_'$i'_R1.fastq.gz' \
$dir_reads'/P004_'$i'_R2.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R1_paired.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R1_unpaired.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R2_paired.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R2_unpaired.fastq.gz' \
ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq3-PE.fa:5:10:20;
done

SAMPLES=(WA05 WB03 WB04 WC02)

for i in ${SAMPLES[@]}; do
java -jar $TRIMMOMATIC/trimmomatic-0.38.jar PE $dir_reads'/P005_'$i'_R1.fastq.gz' \
$dir_reads'/P005_'$i'_R2.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R1_paired.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R1_unpaired.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R2_paired.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R2_unpaired.fastq.gz' \
ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq3-PE.fa:5:10:20;
done

SAMPLES=(WF03 WG03)

for i in ${SAMPLES[@]}; do
java -jar $TRIMMOMATIC/trimmomatic-0.38.jar PE $dir_reads'/P009_'$i'_R1.fastq.gz' \
$dir_reads'/P009_'$i'_R2.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R1_paired.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R1_unpaired.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R2_paired.fastq.gz' \
$dir_trimmed'/grandis_'$i'_R2_unpaired.fastq.gz' \
ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq3-PE.fa:5:10:20;
done

mv *_unpaired.fastq.gz unpaired_reads
```
## Trimming all other *Primula* samples.
Run the following script as follows. It will run the script in pareallel for each sequencing plate.
```sh
for x in {1..9}; do sbatch trimmomatic.sh $x; done
```
Script within trimmomatic.sh
```sh
#!/bin/bash
#SBATCH --job-name=EMC_trimm
#SBATCH -o ./logs/trimm.%j.out
#SBATCH -e ./logs/trimm.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=3800
#SBATCH --nice=1000
#SBATCH --partition=main

TRIMMOMATIC=/home/ubuntu/executables/Trimmomatic-0.38
dir_reads=/home/ubuntu/raw_data
dir_trimmed=/home/ubuntu/trimmed_reads/primula

NUM=${1?Error: no chrom number given}

SAMPLES=$(grep Plate$NUM $dir_trimmed/list_samples_used)

for i in ${SAMPLES[@]}; do
java -jar $TRIMMOMATIC/trimmomatic-0.38.jar PE $dir_reads/$i'_R1.fastq.gz' $dir_reads/$i'_R2.fastq.gz' \
$dir_trimmed/$i'_R1_paired.fastq.gz' \
$dir_trimmed/$i'_R1_unpaired.fastq.gz' \
$dir_trimmed/$i'_R2_paired.fastq.gz' \
$dir_trimmed/$i'_R2_unpaired.fastq.gz' \
ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq3-PE.fa:5:10:20;
done
```
