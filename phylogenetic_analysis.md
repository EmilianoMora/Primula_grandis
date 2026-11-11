# Phylogentic analysis of collinear single copy orthologs
We included sequences from 80 individuals belonging to the six additional species of _Primula_ sect. _Primula_ (_P. elatior_, _P. juliae_, _P. megaseifolia_, _P. renifolia_, _P. veris_, and _P. vulgaris_) and seven outgroup species from five different sections within _Primula_ (_Primula simensis_, _Primula mistassinica_, _Primula souliei_, _Primula nutans_, _Primula boveana_, _Primula verticillata_, and _Dionysia revoluta_).

```sh
#!/bin/bash
#SBATCH --job-name=EM_BCFtools
#SBATCH --error ./logs/BCFtools_merge.%j.err
#SBATCH --out ./logs/BCFtools_merge.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=900:00:00
#SBATCH --mem-per-cpu=3700
#SBATCH --nice=100
#SBATCH --partition=main
#SBATCH --exclude=compute001,med128compute001,med128compute002

BCFTOOLS=/home/ubuntu/executables/bcftools-1.8/bcftools
TABIX=/home/ubuntu/executables/htslib/tabix
VCFTOOLS=/home/ubuntu/executables/vcftools/src/cpp/vcftools
BGZIP=/home/ubuntu/executablesprograms/htslib/bgzip
BEDTOOLS=/home/ubuntu/executablesbedtools2/bin/bedtools
V2P=/home/ubuntu/executables/vcf2phylip.py
PHYX=/home/ubuntu/executables/progs/phyx/src/
IQ=/home/ubuntu/executables/iqtree-2.1.2-Linux/bin/iqtree2
SNPSITES=/home/ubuntu/executables/progs/miniconda3/pkgs/snp-sites-2.5.1-hed695b0_0/bin/snp-sites

# 1) Get the genes we need from the GFF file and create BED file

GFF_FILE=/home/ubuntu/p_grandis_ref_genome/braker_rename.gff
COLLINEAR_BLOCKS=$(ls /home/ubuntu/p_grandis_ref_genome/mcscan_stuff/collinear_blocks/sub_A | grep 'anchor_')

for j in ${COLLINEAR_BLOCKS[@]}; do
LIST_GENES=$(cat /home/ubuntu/p_grandis_ref_genome/mcscan_stuff/collinear_blocks/sub_A/$j | grep -oE 'g[0-9]*.t1')
for i in ${LIST_GENES[@]}; do
grep $i $GFF_FILE | grep CDS | grep + | awk '{print $1,$4,$5}' OFS='\t' >> Coll_Ortho_A.pos.$j.temp.bed
grep $i $GFF_FILE | grep CDS | grep - | awk '{print $1,$4,$5}' OFS='\t' >> Coll_Ortho_A.neg.$j.temp.bed; done
cat Coll_Ortho_A.pos.$j.temp.bed Coll_Ortho_A.neg.$j.temp.bed | sort --version-sort -u -b -k1,1 -k2,2 > Coll_Ortho_A.$j.bed
rm Coll_Ortho_A.pos.$j.temp.bed Coll_Ortho_A.neg.$j.temp.bed; done

mv *bed BED_FILES

# 2) Get the genes we need from the VCF file

VCF_FILE=/home/ubuntu/4_variant_calling/other_primula/subgenome_A/filtering/merge_grandis/merged_allpops_OGs_NoSingletonsDoubletons_variable_biallelic.a.vcf.gz
COLL_ORTHO_BEDS=$(ls BED_FILES | grep bed)

for i in ${COLL_ORTHO_BEDS[@]}; do
$BEDTOOLS intersect -a $VCF_FILE -b ./BED_FILES/$i -header | $BGZIP -c > $i.vcf.gz
$VCFTOOLS --gzvcf $i.vcf.gz --remove eliminate_bad_OGs.txt --max-missing 0.80 --recode --stdout | $BGZIP -c > $i.filtered.vcf.gz
rm $i.vcf.gz
$V2P -i $i.filtered.vcf.gz -f -p #Convert into FASTA file
sed -i '/^[^>]/s/N/\?/g' $i.filtered.min4.fasta #V2P codes heterozygous sites as IUPAC codes. Which are the letters below. I use sed to change those from ambiguity codes to N's remove ambiguous and replace with N
sed '/^[^>]/s/[M|R|W|S|Y|K|V|H|D|B]/N/g' $i.filtered.min4.fasta > $i.filtered.min4.no_ambig.fasta
$SNPSITES $i.filtered.min4.no_ambig.fasta -o $i.filtered.min4.no_ambig.inform.fasta #Now use SNPSITES to remove no informative sites
rm $i.filtered.min4.fasta $i.filtered.min4.no_ambig.fasta; done

mv *vcf.gz VCF_FILES
mv *fasta FASTA_FILES

# 3) IQ tree
#this IQ tree command line will do -B #ultrafast bootstraps, -T threads is set to auto, -m MFP means it will find the best model, -st is telling it we have DNA

$IQ -S ./FASTA_FILES --prefix loci -B 1000 -T AUTO -m MFP+ASC -st DNA -redo \
-o Psiamensis_WA09.SubGenA.bam,Psouliei_WB09.SubGenA.bam,Pbouveana_WC09.SubGenA.bam,Pbouveana_WD09.SubGenA.bam,\
Dionysia_WE08.SubGenA.bam,Pmistassinica_WE09.SubGenA.bam,Pverticilliata_WF08.SubGenA.bam,Pnutans_WF09.SubGenA.bam

$IQ -t /home/ubuntu/emiliano/4th_chp/7_phylogeny/Collinear_Orthologs/sub_A/sing_doubletons/CollOrtho_subgenome_A_OGs_2.treefile --gcf loci.treefile \
-s /home/ubuntu/emiliano/4th_chp/7_phylogeny/Collinear_Orthologs/sub_A/sing_doubletons/Coll_Ortho_subgenome_A_allprim_OGs_2.min4.noambig.inform.fasta --scf 100 --prefix concord_A

mv loci* IQTREE_FILES
mv concord* IQTREE_FILES
```
