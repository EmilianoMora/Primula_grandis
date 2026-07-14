# Varisant calling suing BCFtools
## First *Primula grandis* samples
Run with the following command to run each chromosome in parallel:
```sh
for x in {1..11}; do sbatch BCFtools_calling.sh $x; done
```
Content of the *BCFtools_calling.sh* file
```sh
#!/bin/bash
#SBATCH --job-name=EMC_bcftools
#SBATCH -o ./logs/bcfCall.%j.out
#SBATCH -e ./logs/bcfCall.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=3800
#SBATCH --nice=1000
#SBATCH --partition=main

##set variables
$BCFTOOLS=/home/ubuntu/executables/bcftools-1.8/bcftools

REF=/home/ubuntu/p_grandis_ref_genome/out_JBAT_review4.FINAL_rename.fasta

DIR=/home/ubuntu/variant_calling/p_grandis
BL=$DIR/grandis_bamlist.txt

NUM=${1?Error: no chrom number given}

###Variant calling subgenome A

echo $NUM
$PROGRAMS/bcftools-1.8/bcftools mpileup -r scaf_${NUM}a -b $BL --threads $SLURM_CPUS_PER_TASK -a AD,DP,SP -f $REF -Ou | $BCFTOOLS call -m -f GQ,GP \
--threads $SLURM_CPUS_PER_TASK -Oz -o $DIR/pgrandis_${NUM}.a.vcf.gz

##then need to index
echo "Finished pileup and calling. Now index"
$BCFTOOLS index -f $DIR/pgrandis_${NUM}.a.vcf.gz

###Variant calling subgenome B

echo $NUM
$PROGRAMS/bcftools-1.8/bcftools mpileup -r scaf_${NUM}b -b $BL --threads $SLURM_CPUS_PER_TASK -a AD,DP,SP -f $REF -Ou | $BCFTOOLS call -m -f GQ,GP \
--threads $SLURM_CPUS_PER_TASK -Oz -o $DIR/pgrandis_${NUM}.b.vcf.gz

##then need to index
echo "Finished pileup and calling. Now index"
$BCFTOOLS index -f $DIR/pgrandis_${NUM}.b.vcf.gz
```
## For other *Primula*
Run with the following command to run each chromosome in parallel:
```sh
for x in {1..11}; do sbatch BCFtools_calling.sh $x; done
```
Content of the *BCFtools_calling.sh* file
```sh
#!/bin/bash
#SBATCH --job-name=EMC_bcftools
#SBATCH -o ./logs/bcfCall.%j.out
#SBATCH -e ./logs/bcfCall.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=7700
#SBATCH --nice=1000
#SBATCH --partition=main

NUM=${1?Error: no chrom number given}

##set variables
$BCFTOOLS=/home/ubuntu/executables/bcftools-1.8/bcftools

REF=/home/ubuntu/p_grandis_ref_genome/out_JBAT_review4.FINAL_rename.fasta

DIR=$DIR/variant_calling/other_primula

#Variant calling subgenome A

BL=$DIR/bamlist.txt

echo $NUM

$PROGRAMS/bcftools-1.8/bcftools mpileup -r scaf_${NUM}a -b $BL --threads $SLURM_CPUS_PER_TASK -a AD,DP,SP -f $REF -Ou | $BCFTOOLS call -m -f GQ,GP \
--threads $SLURM_CPUS_PER_TASK -Oz -o $DIR/subgenome_A/primula_${NUM}.a.vcf.gz

$BCFTOOLS reheader --threads $SLURM_CPUS_PER_TASK --samples $DIR/rename_samples_A.txt -o $DIR/primula_${NUM}_rehead.a.vcf.gz $DIR/subgenome_A/primula_${NUM}.a.vcf.gz

rm $DIR/subgenome_A/primula_${NUM}.a.vcf.gz

##then need to index
echo "Finished pileup and calling. Now index"
$BCFTOOLS index $DIR2/subgenome_A/primula_${NUM}_rehead.a.vcf.gz

#Variant calling subgenome B

echo $NUM

$PROGRAMS/bcftools-1.8/bcftools mpileup -r scaf_${NUM}a -b $BL --threads $SLURM_CPUS_PER_TASK -a AD,DP,SP -f $REF -Ou | $BCFTOOLS call -m -f GQ,GP \
--threads $SLURM_CPUS_PER_TASK -Oz -o $DIR/subgenome_B/primula_${NUM}.b.vcf.gz

$BCFTOOLS reheader --threads $SLURM_CPUS_PER_TASK --samples $DIR/rename_samples_B.txt -o $DIR/primula_${NUM}_rehead.b.vcf.gz $DIR/subgenome_B/primula_${NUM}.b.vcf.gz

rm $DIR/subgenome_B/primula_${NUM}.b.vcf.gz

##then need to index
echo "Finished pileup and calling. Now index"
$BCFTOOLS index $DIR/subgenome_B/primula_${NUM}_rehead.b.vcf.gz
```
