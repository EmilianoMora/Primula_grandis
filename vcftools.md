# Code to perform filtering of VCF file
Run with the following command to run all 11 chromosomes in pareallel:
```sh
for x in {1..11}; do sbatch FILTERING.sh $x; done
```
Content of the *FILTERING.sh* file:
```sh
#!/bin/bash
#SBATCH --job-name=EM_VCFtools
#SBATCH --error ./logs/filter_variants.%j.err
#SBATCH --out ./logs/filter_variants.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=3700

BCFTOOLS=/home/ubuntu/executables/bcftools-1.8/bcftools
VCFTOOLS=/home/ubuntu/executables/vcftools/src/cpp/vcftools
BGZIP=/home/ubuntu/executables/htslib/bgzip
TABIX=/home/ubuntu/executables/htslib/tabix

NUM=${1?Error: no chrom number given}

echo $NUM

#Filtering chromosome A

INPUT=../pgrandis_$NUM.a.vcf.gz
OUTPUT_1=pgrandis_$NUM.a.no_indel.vcf.gz
OUTPUT_2=pgrandis_$NUM.a.no_indel.filtered.vcf.gz
OUTPUT_3=pgrandis_$NUM.a.no_indel.filtered.AB.vcf.gz
OUTPUT_4=pgrandis_$NUM.a.no_indel.filtered.AB.rh.vcf.gz

# set filters parameters
#MISS=0.8
QUAL=30
MIN_DEPTH=5
MAX_DEPTH=60
MAX_ALLELE=2

#Basic filtering

$BCFTOOLS view -e 'IMF > 0.0001' $INPUT --threads $SLURM_CPUS_PER_TASK -Oz -o $OUTPUT_1

$VCFTOOLS --gzvcf $OUTPUT_1 --remove-indels --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --max-alleles $MAX_ALLELE --recode --recode-INFO-all --stdout | $BGZIP -c > $OUTPUT_2

rm $OUTPUT_1

#Filter allelic imbalance

$BCFTOOLS filter -e 'GT="het" & (FORMAT/AD[*:1])/((FORMAT/AD[*:0]) + (FORMAT/AD[*:1]))  <= 0.30 | GT="het" & (FORMAT/AD[*:1])/((FORMAT/AD[*:0]) + (FORMAT/AD[*:1])) >= 0.70' \
$OUTPUT_2 --threads $SLURM_CPUS_PER_TASK -Oz -o $OUTPUT_3

rm $OUTPUT_2

#Change sample names (currently they all start with '/home/ubuntu/emiliano/4th_chp/3_mapped_reads/p_grandis/*.bam')

$BCFTOOLS view --header $OUTPUT_3 > header.$NUM.txt
sed -i 's|/home/ubuntu/emiliano/4th_chp/3_mapped_reads/p_grandis/||g' header.$NUM.txt

$BCFTOOLS reheader -h header.$NUM.txt $OUTPUT_3 -o $OUTPUT_4

rm $OUTPUT_3 header.$NUM.txt
$TABIX -p vcf $OUTPUT_4

$VCFTOOLS --gzvcf $OUTPUT_4 --het --out het/$NUM.a

#######################################################################################################
#Filtering chromosome B
#######################################################################################################

INPUT=../pgrandis_$NUM.b.vcf.gz
OUTPUT_1=pgrandis_$NUM.b.no_indel.vcf.gz
OUTPUT_2=pgrandis_$NUM.b.no_indel.filtered.vcf.gz
OUTPUT_3=pgrandis_$NUM.b.no_indel.filtered.AB.vcf.gz
OUTPUT_4=pgrandis_$NUM.b.no_indel.filtered.AB.rh.vcf.gz

# set filters parameters
MISS=0.8
QUAL=30
MIN_DEPTH=5
MAX_DEPTH=60
MAX_ALLELE=2

#Basic filtering

$BCFTOOLS view -e 'IMF > 0.0001' $INPUT --threads $SLURM_CPUS_PER_TASK -Oz -o $OUTPUT_1

$VCFTOOLS --gzvcf $OUTPUT_1 --remove-indels --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --max-alleles $MAX_ALLELE --recode --recode-INFO-all --stdout | $BGZIP -c > $OUTPUT_2

rm $OUTPUT_1

#Filter allelic imbalance

$BCFTOOLS filter -e 'GT="het" & (FORMAT/AD[*:1])/((FORMAT/AD[*:0]) + (FORMAT/AD[*:1]))  <= 0.30 | GT="het" & (FORMAT/AD[*:1])/((FORMAT/AD[*:0]) + (FORMAT/AD[*:1])) >= 0.70' \
$OUTPUT_2 --threads $SLURM_CPUS_PER_TASK -Oz -o $OUTPUT_3

rm $OUTPUT_2

#Change sample names (currently they all start with '/home/ubuntu/emiliano/4th_chp/3_mapped_reads/p_grandis/*.bam')

$BCFTOOLS view --header $OUTPUT_3 > header.$NUM.txt
sed -i 's|/home/ubuntu/emiliano/4th_chp/3_mapped_reads/p_grandis/||g' header.$NUM.txt

$BCFTOOLS reheader -h header.$NUM.txt $OUTPUT_3 -o $OUTPUT_4

rm $OUTPUT_3 header.$NUM.txt
$TABIX -p vcf $OUTPUT_4

#$VCFTOOLS --gzvcf $OUTPUT_4 --het --out het/$NUM.b
```
# Concatenate all filtered VCF file sinto one
```sh
#!/bin/bash
#SBATCH --job-name=EM_gatk_filter
#SBATCH --error ./logs/concatenate_vcf.%j.err
#SBATCH --out ./logs/concatenate_vcf.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=3800
#SBATCH --nice=0
#SBATCH --partition=main

BCFTOOLS=/home/ubuntu/executables/bcftools-1.8/bcftools
VCFTOOLS=/home/ubuntu/executables/vcftools/src/cpp/vcftools
BGZIP=/home/ubuntu/executables/htslib/bgzip
TABIX=/home/ubuntu/executables/htslib/tabix

#Concatenate all the chromosomes into one VCF file.

$BCFTOOLS concat -f list_chr.txt -Oz --threads $SLURM_CPUS_PER_TASK > pgrandis_concat.vcf.gz
$TABIX -p vcf pgrandis_concat.vcf.gz

#Concatenate per subgenome. This one ends up with a VCF file of all samples for each chromosome

$BCFTOOLS concat -f list_chr_subgen_A.txt -Oz --threads $SLURM_CPUS_PER_TASK > pgrandis_concat_subgen_A.vcf.gz
$TABIX -p vcf pgrandis_concat_subgen_A.vcf.gz

$BCFTOOLS concat -f list_chr_subgen_B.txt -Oz --threads $SLURM_CPUS_PER_TASK > pgrandis_concat_subgen_B.vcf.gz
$TABIX -p vcf pgrandis_concat_subgen_B.vcf.gz
```
