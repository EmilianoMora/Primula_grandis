# Estimation genomewide divergence 

To estimate approximate time of divergence of both subgenomes (A and B) we estimated divergence at synonymous sites between collinear orthologs of the two subgenomes using (ParaAT2.0)[https://github.com/jdebarry/paraat/tree/master]. The collinear orthologs between the subgenomes were identified using [MCscan](https://github.com/tanghaibao/jcvi/wiki/Mcscan-(python-version). For this analyses both the DNA and protein fasta files are needed.

```sh
#!/bin/bash
#SBATCH --job-name=EMC_kaks
#SBATCH --output=./logs/kaks.%j.out
#SBATCH --error ./logs/kaks.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=3700
#SBATCH --partition=main
#SBATCH --nice=1000

DIR=/home/ubuntu/4th_chp/9_dNdS/pgraA_pgraB
rm $DIR/pgra.longCDS.Collinear_Orthologs.fasta $DIR/pgra.longProt.Collinear_Orthologs.fasta
rm -r $DIR/pgra_kaks

#Step 1) Create necessary input files for ParaAT

GENE=$(awk '1' RS='[[:space:]]+' $DIR/files/pgraA.pgraB.anchors)

for i in ${GENE[@]}; do
grep -A1 $i $DIR/files/pgra.longCDS_lin.fasta >> $DIR/pgra.longCDS.Collinear_Orthologs.fasta
sed -n '/'$i'/,/>/p' $DIR/files/pgra.longProt.fasta | head -n -1 >> $DIR/pgra.longProt.Collinear_Orthologs.fasta;done

#Step 2) run ParaAT
##  SET THE VARIABLES BELOW  ######################
GRANDIS=$DIR/files/pgraA.pgraB.anchors
CDS=$DIR/pgra.longCDS.Collinear_Orthologs.fasta
PROT=$DIR/pgra.longProt.Collinear_Orthologs.fasta
###################################################

export PATH=/home/ubuntu/programs/ParaAT2.0:$PATH; \
export PATH=/home/ubuntu/programs/KaKs_Calculator2.0/bin/Linux:$PATH; \
export PATH=/home/ubuntu/programs:$PATH; \
~/programs/ParaAT2.0/ParaAT.pl \
-h $GRANDIS \
-n $CDS \
-a $PROT \
-o $DIR/pgrandis_kaks \
-p processors \
-m muscle \
-f axt \
-k \
-v

#Step 3) make summary file of kaks
FILE=$(ls $DIR/pgrandis_kaks | grep cds_aln.axt.kaks)
for i in ${FILE[@]}; do sed -n '2p' $DIR/pgrandis_kaks/$i >> $DIR/pgrandis_kaks/summary_kaks_grandisAB;done
```
