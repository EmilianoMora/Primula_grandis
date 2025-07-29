# TEST
## Subgenome characterization of *Primula grandis* with SubPhaser
If P. grandis is allotetraploid, as proposed by previous phylogenetic analyses (Schmidt-Lebuhn et al., 2012; Stubbs et al., 2023), it should contain two distinct subgenomes. We tested this prediction by using subgenome-specific k-mers-analyses of LTRs implemented in SubPhaser v.1.2.6 (Jia et al., 2022).
```sh
conda activate SubPhaser
subphaser -i out_JBAT_review4.FINAL_22scaf_rename.fa.gz -c sg.config -p 32 -max_memory 7.8G
conda deactivate
```
The content of sg.config is the following:
```
scaf_1a scaf_1b
scaf_2a scaf_2b
scaf_3a scaf_3b
scaf_4a scaf_4b
scaf_5a scaf_5b
scaf_6a scaf_6b
scaf_7a scaf_7b
scaf_8a scaf_8b
scaf_9a scaf_9b
scaf_10a scaf_10b
scaf_11a scaf_11b
```
## Subgenome dominance in *P. grandis*
### Mapping RNA sequencing reads to reference genome using Bowtie
To test for subgenome dominance, we checked for Homoeologous Expression Bias (HEB) towards one of the subgenomes of P. grandis as follows. First, we aligned transcriptomic data from leaves, flowers and floral buds to the reference genome using Bowtie2 v2.3.4.3 (Langmead & Salzberg, 2012) with the -very-sensitive-local option (-D 20 -R 3 -N 0 -L 20 -i S, 1, 0.50).
```sh
#!/bin/bash
#SBATCH --job-name=bowtie
#SBATCH --output=./logs/bowtie.%A_%a.out
#SBATCH --error=./logs/bowtie.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=500:00:00
#SBATCH --mem-per-cpu=3800
#SBATCH --nice=1000
#SBATCH --partition=main
#SBATCH --array=1-3%3

# ==== CONFIGURATION ====
# Reference index basename (built with bowtie2-build)
INDEX="out_JBAT_review4.FINAL_rename"

# Path to trimmed paired-end FASTQ files
RNA_DATA="/home/ubuntu/my_data-2/analysis/map/grand/RNA_seq/trimmomatic"

# Tools
BOWTIE="/home/ubuntu/emiliano/executables/bowtie2-2.3.4.3-linux-x86_64/bowtie2"
SAMTOOLS="/home/ubuntu/emiliano/executables/samtools-1.10/samtools"

# ==== GET INPUT FILES ====
FORWARD=$(ls "$RNA_DATA"/*.qc.R1.paired.fastq | sed -n "${SLURM_ARRAY_TASK_ID}p")
REVERSE="${FORWARD/R1/R2}"
NAME=$(basename "$FORWARD" .qc.R1.paired.fastq)
OUTPUT="./${NAME}.bam"

# ==== ALIGNMENT ====
echo "Starting alignment for $NAME"
"$BOWTIE" -p "$SLURM_CPUS_PER_TASK" -q --very-sensitive-local -x "$INDEX" -1 "$FORWARD" -2 "$REVERSE" -S "${NAME}.sam"

# ==== CONVERT SAM TO BAM, SORT ====
echo "Converting and sorting alignment for $NAME"
"$SAMTOOLS" view -@ "$SLURM_CPUS_PER_TASK" -b "${NAME}.sam" | \
    "$SAMTOOLS" sort -@ "$SLURM_CPUS_PER_TASK" -T "${NAME}.tmp" -o "$OUTPUT"

# ==== INDEX BAM FILE ====
echo "Indexing BAM file for $NAME"
"$SAMTOOLS" index "$OUTPUT"

# ==== CLEANUP ====
rm -f "${NAME}.sam"
echo "Finished processing $NAME at $(date)"
```
### Estimating Reads Per Kilobase per Million (RPKM)
Secondly, the number of aligned reads mapping to each gene was counted using htseq-count v2.0.5 (Putri et al., 2022) to estimate Reads Per Kilobase per Million (RPKM). Thirdly, given that we cannot use the likelihood ratio test proposed by Smith et al. (2019) due to the lack of biological and tissue replicates, we established -3 ≤ log2 HEB ≥ 3 (i.e., 8-fold HEB bias towards subgenome-A or -B, respectively) as the threshold of HEB significance.
