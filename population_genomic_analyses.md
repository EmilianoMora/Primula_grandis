# Estimates of genetic diversity
To test the expectations of genome-wide reduced genetic diversity associated with both polyploid-speciation and shift to self-compatibility, we compared genome-wide estimates of genetic diversity between both P. grandis subgenomes and the six diploid species in the Caucasus region (P. elatior, P. juliae, P. megaseifolia, P. renifolia, P. veris, and P. vulgaris). To do this, we estimated the number of segregating sites (S), and nucleotide diversity (π) at synonymous and non-synonymous sites using the python script ‘popgenWindows.py’ (https://github.com/simonhmartin/genomics_general) in nonoverlapping 50-kb windows across the genome. Moreover, to determine whether P. grandis subgenomes have genome-wide patterns consistent with recent changes in effective population size, we compared the genome-wide distributions of Tajima’s D (D) of P. grandis subgenomes-A and -B and all its relatives in Primula sect. Primula. Tajima’s D was calculated using nonoverlapping 50-kb windows using the ‘popgenWindows.py’ python script as above.
## Estimates of *P. grandis* subgenomes
```sh
# Calculate individual heterozygosity only on chr 01, using edited windows that exclude the S locus.
# To calculate individual heterozygosity on chr 02-11, use run_sm_popgenwindows_indHet.sh.

SM=/home/ubuntu/emiliano/BACKUP_STUBBS/scripts_rebecca_2/simonmartin_fst/genomics_general
DIR=/home/ubuntu/emiliano/4th_chp/4_variant_calling/p_grandis/annotate_sites
FILES=(0fold_pgrandis_A 0fold_pgrandis_B 4fold_pgrandis_A 4fold_pgrandis_B)
for i in ${FILES[@]}; do
NAME=$i
VCF=$i.vcf.gz
GENO=$i.geno.gz
POP=popFile.txt
echo "### Create .geno file for $i"
python3 $SM/VCF_processing/parseVCF.py -i $DIR/$VCF -o $GENO
echo "### Calculate Fst and Dxy for $i"
python3 $SM/popgenWindows.py -g $GENO -o $NAME.popgenWindows.pi.grandis.csv.gz -f phased -w 50000 -m 100 \
-T $SLURM_CPUS_PER_TASK --addWindowID --writeFailedWindows --verbose -p A --popsFile $POP
python3 $SM/popgenWindows.py -g $GENO -o $NAME.popgenWindows.Tajima.csv.gz -f phased -w 50000 -m 100 \
-T $SLURM_CPUS_PER_TASK --addWindowID --writeFailedWindows --verbose -p A --analysis popFreq --popsFile $POP; done

#echo "### Calculate heterozygosity per each individual"
python3 $SM/popgenWindows.py -g $GENO \
-o $OUTFOLDER/$NAME.popgenWindows_$WSIZE\k.indHet.csv.gz \
-f phased -w $WSIZE\000 -m 100 -T $SLURM_CPUS_PER_TASK --addWindowID --writeFailedWindows --verbose -p pin -p thrum \
--analysis indHet \
--popsFile $SETS

echo "done?"
```
Content of popFile.txt
```grandis_WA05.bam	A
grandis_WB03.bam	A
grandis_WB04.bam	A
grandis_WB10.bam	A
grandis_WC02.bam	A
grandis_WC09.bam	A
grandis_WD08.bam	A
grandis_WE08.bam	A
grandis_WF03.bam	A
grandis_WG03.bam	A
```
## Estimation of genetic diversity in other species within *Primula* sect. *Primula*
```sh
# Calculate individual heterozygosity only on chr 01, using edited windows that exclude the S locus.
# To calculate individual heterozygosity on chr 02-11, use run_sm_popgenwindows_indHet.sh.

#SM=/home/ubuntu/rebecca-3/simonmartin_fst/genomics_general #OLD DIRECTORY
SM=/home/ubuntu/emiliano/BACKUP_STUBBS/scripts_rebecca_2/simonmartin_fst/genomics_general
DIR=/home/ubuntu/emiliano/4th_chp/4_variant_calling/other_primula/subgenome_B/filtering/annotate
FILES=(0fold 4fold)

for i in ${FILES[@]}; do
NAME=$i
VCF=''$i'_merged_sect_primula_short.b.vcf.gz'
GENO=$i.geno.gz
POP=sect_primula_popgen.txt
echo "### Create .geno file for $i"
python3 $SM/VCF_processing/parseVCF.py -i $DIR/$VCF -o $GENO
echo "### Calculate Fst and Dxy for $i"
python3 $SM/popgenWindows.py -g $GENO -o $NAME.popgenWindows.pi.dxy.Fst.sect_primula.csv.gz -f phased -w 50000 -m 100 \
-T $SLURM_CPUS_PER_TASK --addWindowID --writeFailedWindows --verbose -p pelatior -p pmega -p preni -p pveris -p pvulg -p pgrand_B --popsFile $POP
echo "### Calculate Tajima's D for $i"
python3 $SM/popgenWindows.py -g $GENO -o $NAME.popgenWindows.Tajima.sect_primula.csv.gz -f phased -w 50000 -m 100 \
-T $SLURM_CPUS_PER_TASK --addWindowID --writeFailedWindows --verbose -p pelatior -p pmega -p preni -p pveris -p pvulg -p pgrand_B --analysis popFreq --popsFile $POP; done
```
Content of sect_primula_popgen.txt
```
Pela_palla_WD02.SubGenB.bam	pelatior
Pela_palla_WC03.SubGenB.bam	pelatior
Pela_palla_WH04.SubGenB.bam	pelatior
Pela_palla_WB02.SubGenB.bam	pelatior
Pela_palla_WA03.SubGenB.bam	pelatior
Pela_palla_WB01.SubGenB.bam	pelatior
Pver_macro_WF06.SubGenB.bam	pveris
Pver_macro_WE07.SubGenB.bam	pveris
Pver_macro_WA05.SubGenB.bam	pveris
Pver_macro_WH05.SubGenB.bam	pveris
Pver_macro_WG06.SubGenB.bam	pveris
Pver_macro_WF07.SubGenB.bam	pveris
Pver_macro_WH04.SubGenB.bam	pveris
Pver_macro_WG05.SubGenB.bam	pveris
Pmega_WC06.SubGenB.bam	pmega
Pmega_WD03.SubGenB.bam	pmega
Pmega_WD12.SubGenB.bam	pmega
Pmega_WC07.SubGenB.bam	pmega
Pmega_WE11.SubGenB.bam	pmega
Pmega_WD02.SubGenB.bam	pmega
Preni_WH01.SubGenB.bam	preni
Preni_WH02.SubGenB.bam	preni
Preni_WH05.SubGenB.bam	preni
Preni_WA05.SubGenB.bam	preni
Preni_WA06.SubGenB.bam	preni
Preni_WA07.SubGenB.bam	preni
Preni_WA08.SubGenB.bam	preni
Preni_WA09.SubGenB.bam	preni
Preni_WA10.SubGenB.bam	preni
Preni_WA11.SubGenB.bam	preni
Pvul_vul_WG02.SubGenB.bam	pvulg
Pvul_vul_WG09.SubGenB.bam	pvulg
Pvul_vul_WG10.SubGenB.bam	pvulg
Pvul_vul_WF11.SubGenB.bam	pvulg
Pvul_vul_WE12.SubGenB.bam	pvulg
Pvul_vul_WE01.SubGenB.bam	pvulg
Pvul_vul_WE02.SubGenB.bam	pvulg
Pvul_vul_WE03.SubGenB.bam	pvulg
Pvul_vul_WE04.SubGenB.bam	pvulg
Pvul_vul_WE05.SubGenB.bam	pvulg
Pvul_vul_WE06.SubGenB.bam	pvulg
```
## Code for R plots
```R
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggridges)
library(data.table)

setwd("~/r_analysis/primula_grandis/databases/popgen")

four_fold_pi_sect_primula <- read.csv("4fold.popgenWindows.pi.dxy.Fst.sect_primula.csv", header = T)
four_fold_pi_grandis_A <- read.csv("4fold_pgrandis_A.popgenWindows.pi.grandis.csv", header = T)
four_fold_pi_grandis_B <- read.csv("4fold_pgrandis_B.popgenWindows.pi.grandis.csv", header = T)

# Nucleotide diversity (pi)

summary(four_fold_pi_sect_primula$pi_pelatior)
summary(four_fold_pi_sect_primula$pi_pmega)
summary(four_fold_pi_sect_primula$pi_preni)
summary(four_fold_pi_sect_primula$pi_pveris)
summary(four_fold_pi_sect_primula$pi_pvulg)
summary(four_fold_pi_sect_primula$pi_pgrand)

summary(four_fold_pi_grandis_A$pi_A)
sd(four_fold_pi_grandis_A$pi_A,na.rm = T)
summary(four_fold_pi_grandis_B$pi_A)
sd(four_fold_pi_grandis_B$pi_A,na.rm = T)

############################
#PLOT NUCLEOTIDE DIVERISITY
############################

four_fold_pi_sect_primula_melted <- melt(four_fold_pi_sect_primula, id='windowID')  %>%  
  filter(variable %in% c('pi_pelatior','pi_pmega','pi_preni','pi_pveris','pi_pvulg'))

colnames(four_fold_pi_grandis_A)[7] <- "pi_pgrand_sub_A" #change name of column 7
four_fold_pi_grandis_A_melted <- melt(four_fold_pi_grandis_A,id='windowID')  %>%  
  filter(variable == 'pi_pgrand_sub_A')

colnames(four_fold_pi_grandis_B)[7] <- "pi_pgrand_sub_B" #change name of column 7
four_fold_pi_grandis_B_melted <- melt(four_fold_pi_grandis_B,id='windowID')  %>%  
  filter(variable == 'pi_pgrand_sub_B')

#merge the three dataframes (sect_primula, grandis A and grandis B)

df_list <- list(four_fold_pi_sect_primula_melted, four_fold_pi_grandis_A_melted, four_fold_pi_grandis_B_melted)
four_fold_pi_sect_primula_melted_2 <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)

colnames(four_fold_pi_sect_primula_melted_2) <- c('windowID','species','pi')
four_fold_pi_sect_primula_melted_2$pi <-as.numeric(four_fold_pi_sect_primula_melted_2$pi)

p <- ggplot(four_fold_pi_sect_primula_melted_2, aes(x=species, y=pi, fill=species)) + 
  geom_violin(trim = T, na.rm = T, scale = "width", colour = "white") + 
  geom_boxplot(width=0.30, fill="white", na.rm = T, outlier.colour = "#565656", outlier.size = 0.6) + 
  ylim(0,0.015) + theme_classic(base_size = 11) +
  labs(title="", x="", y = "Nucleotide Diversity at Synonymous Sites") + 
  scale_fill_manual(values=c("#ff7f00", "#75CCCB", "#FFE733","#4daf4a", "#f781bf", "#C73333","#377eb8"), 
                    labels=c("P. elatior","P. megaseifolia","P. renifolia","P. veris","P. vulgaris","P. grandis_A","P. grandis_B")) +
  scale_x_discrete(labels=c("P. elatior","P. megaseifolia","P. renifolia","P. veris","P. vulgaris","P. grandis_A","P. grandis_B")) +
  theme(axis.text.x = element_text(face = "italic", size=11, color = "black")) +
  theme(axis.text.y = element_text(size=10, color = "black"))

##########################
#piN/piS
##########################

setwd("~/r_analysis/primula_grandis/databases/popgen")

four_fold_pi_sect_primula <- read.csv("4fold.popgenWindows.pi.dxy.Fst.sect_primula.csv", header = T)
four_fold_pi_sect_primula <- four_fold_pi_sect_primula[ -c(13:42)] %>%  #Eliminating unnecessary columns
  rename_with(~paste0(., "_4fold"), grep("^[A-Z]*", names(.))) #add "_4fold" to all column names

zero_fold_pi_sect_primula <- read.csv("0fold.popgenWindows.pi.dxy.Fst.sect_primula.csv", header = T)
zero_fold_pi_sect_primula <- zero_fold_pi_sect_primula[ -c(13:42)] %>%  #Eliminating unnecessary columns
  rename_with(~paste0(., "_0fold"), grep("^[A-Z]*", names(.))) #add "_4fold" to all column names

sect_primula <- cbind(four_fold_pi_sect_primula, zero_fold_pi_sect_primula)

sect_primula$pi_pelatior_4fold <- ifelse(sect_primula$pi_pelatior_4fold<=0.0001,0,sect_primula$pi_pelatior_4fold) # Turn everything that is below 0.0009 into a zero. Low values of piN drags piN/piS upwards.
sect_primula$pi_pelatior_0fold <- ifelse(sect_primula$pi_pelatior_0fold<=0.0005,0,sect_primula$pi_pelatior_0fold) 

sect_primula$pi_preni_4fold <- ifelse(sect_primula$pi_preni_4fold<=0.0001,0,sect_primula$pi_preni_4fold) # Turn everything that is below 0.0009 into a zero. Low values of piN drags piN/piS upwards.
sect_primula$pi_preni_0fold <- ifelse(sect_primula$pi_preni_0fold<=0.0005,0,sect_primula$pi_preni_0fold) 

sect_primula$pi_pmega_4fold <- ifelse(sect_primula$pi_pmega_4fold<=0.0001,0,sect_primula$pi_pmega_4fold) # Turn everything that is below 0.0009 into a zero. Low values of piN drags piN/piS upwards.
sect_primula$pi_pmega_0fold <- ifelse(sect_primula$pi_pmega_0fold<=0.0005,0,sect_primula$pi_pmega_0fold) 

sect_primula$pi_pveris_4fold <- ifelse(sect_primula$pi_pveris_4fold<=0.0001,0,sect_primula$pi_pveris_4fold) # Turn everything that is below 0.0009 into a zero. Low values of piN drags piN/piS upwards.
sect_primula$pi_pveris_0fold <- ifelse(sect_primula$pi_pveris_0fold<=0.0005,0,sect_primula$pi_pveris_0fold) 

sect_primula$pi_pvulg_4fold <- ifelse(sect_primula$pi_pvulg_4fold<=0.0001,0,sect_primula$pi_pvulg_4fold) # Turn everything that is below 0.0009 into a zero. Low values of piN drags piN/piS upwards.
sect_primula$pi_pvulg_0fold <- ifelse(sect_primula$pi_pvulg_0fold<=0.0005,0,sect_primula$pi_pvulg_0fold) 

sect_primula <- sect_primula %>% mutate(piN_piS_pela = sect_primula$pi_pelatior_0fold/sect_primula$pi_pelatior_4fold) %>%
  mutate(piN_piS_pmega = sect_primula$pi_pmega_0fold/sect_primula$pi_pmega_4fold) %>%
  mutate(piN_piS_preni = sect_primula$pi_preni_0fold/sect_primula$pi_preni_4fold) %>%
  mutate(piN_piS_pveris = sect_primula$pi_pveris_0fold/sect_primula$pi_pveris_4fold) %>%
  mutate(piN_piS_pvulg = sect_primula$pi_pvulg_0fold/sect_primula$pi_pvulg_4fold)
#  mutate(piN_piS_pgrand_B = sect_primula$pi_pgrand_B_0fold/sect_primula$pi_pgrand_B_4fold)

sect_primula_reduced <- sect_primula[-c(2:24)]
sect_primula_melted <- reshape2::melt(sect_primula_reduced, id='windowID_4fold')
colnames(sect_primula_melted) <- c('windowID','species','piN_piS')

four_fold_pi_grandis_A <- read.csv("4fold_pgrandis_A.popgenWindows.pi.grandis.csv", header = T)
zero_fold_pi_grandis_A <- read.csv("0fold_pgrandis_A.popgenWindows.pi.grandis.csv", header = T)
grandis_A <- cbind(four_fold_pi_grandis_A, zero_fold_pi_grandis_A)
grandis_A <- grandis_A[-c(2:6,8:13)]
#grandis_A <- grandis_A %>%  filter(pi_A > 0, pi_A.1 > 0)
grandis_A <- grandis_A %>% mutate(piN_piS_pgrand_A = grandis_A$pi_A.1/grandis_A$pi_A)

four_fold_pi_grandis_B <- read.csv("4fold_pgrandis_B.popgenWindows.pi.grandis.csv", header = T)
zero_fold_pi_grandis_B <- read.csv("0fold_pgrandis_B.popgenWindows.pi.grandis.csv", header = T)
grandis_B <- cbind(four_fold_pi_grandis_B,zero_fold_pi_grandis_B)
grandis_B <- grandis_B[-c(2:6,8:13)]
#grandis_B <- grandis_B %>%  filter(pi_A > 0, pi_A.1 > 0)
grandis_B <- grandis_B %>% mutate(piN_piS_pgrand_B = grandis_B$pi_A.1/grandis_B$pi_A)

grandis_A_B <- merge(data.frame(grandis_A, row.names=NULL), data.frame(grandis_B, row.names=NULL), by = 0, all = TRUE)[-1] #merge two dataframes even when they don't have same number of rows
grandis_A_B <- grandis_A_B[-c(2,3,5:7)]
###problem is here
grandis_A_B_melted <- reshape2::melt(grandis_A_B, id='windowID.x')
colnames(grandis_A_B_melted) <- c('windowID','species','piN_piS')

all_together <- rbind(sect_primula_melted, grandis_A_B_melted)
is.na(all_together) <- sapply(all_together, is.infinite)

all_together %>% group_by(species) %>% summarize(mean=mean(piN_piS, na.rm = TRUE), 
                                                 sd=sd(piN_piS, na.rm = TRUE), 
                                                 median=median(piN_piS, na.rm = TRUE))

ggplot(all_together, aes(x=species, y=piN_piS, fill=species)) + 
  geom_violin(trim = T, na.rm = T, scale = "width", colour = "white") + 
  geom_boxplot(width=0.30, fill="white", na.rm = T, outlier.colour = "#565656", outlier.size = 0.6) + 
  ylim(0,5) + theme_classic(base_size = 11) +
  labs(title="", x="", y = "pi-Non-Synonymous / pi-Synonymous") + 
  scale_fill_manual(values=c("#ff7f00", "#75CCCB", "#FFE733","#4daf4a", "#f781bf", "#C73333","#377eb8"), 
                    labels=c("P. elatior","P. megaseifolia","P. renifolia","P. veris","P. vulgaris","P. grandis_A","P. grandis_B")) +
  scale_x_discrete(labels=c("P. elatior","P. megaseifolia","P. renifolia","P. veris","P. vulgaris","P. grandis_A","P. grandis_B")) +
  theme(axis.text.x = element_text(face = "italic", size=11, color = "black")) +
  theme(axis.text.y = element_text(size=10, color = "black"))

######

four_fold_pi_sect_primula <- read.csv("4fold.popgenWindows.pi.dxy.Fst.sect_primula.csv", header = T)
four_fold_pi_grandis_A <- read.csv("4fold_pgrandis_A.popgenWindows.pi.grandis.csv", header = T)
four_fold_pi_grandis_B <- read.csv("4fold_pgrandis_B.popgenWindows.pi.grandis.csv", header = T)

zero_fold_pi_sect_primula <- read.csv("0fold.popgenWindows.pi.dxy.Fst.sect_primula.csv", header = T)
zero_fold_pi_grandis_A <- read.csv("0fold_pgrandis_A.popgenWindows.pi.grandis.csv", header = T)
zero_fold_pi_grandis_B <- read.csv("0fold_pgrandis_B.popgenWindows.pi.grandis.csv", header = T)

four_fold_pi_sect_primula$pi_pelatior <- ifelse(four_fold_pi_sect_primula$pi_pelatior<=0.0005,0,four_fold_pi_sect_primula$pi_pelatior) #Change nucleotide diveristy values below 0.0001 as '0'
four_fold_pi_sect_primula$pi_pmega <- ifelse(four_fold_pi_sect_primula$pi_pmega<=0.0005,0,four_fold_pi_sect_primula$pi_pmega)
four_fold_pi_sect_primula$pi_preni <- ifelse(four_fold_pi_sect_primula$pi_preni<=0.0005,0,four_fold_pi_sect_primula$pi_preni)
four_fold_pi_sect_primula$pi_pveris <- ifelse(four_fold_pi_sect_primula$pi_pveris<=0.0005,0,four_fold_pi_sect_primula$pi_pveris)
four_fold_pi_sect_primula$pi_pvulg <- ifelse(four_fold_pi_sect_primula$pi_pvulg<=0.0005,0,four_fold_pi_sect_primula$pi_pvulg)
four_fold_pi_sect_primula$pi_pgrand_B <- ifelse(four_fold_pi_sect_primula$pi_pgrand_B<=0.0005,0,four_fold_pi_sect_primula$pi_pgrand_B)

zero_fold_pi_sect_primula$pi_pelatior <- ifelse(zero_fold_pi_sect_primula$pi_pelatior<=0.0005,0,zero_fold_pi_sect_primula$pi_pelatior) #Change nucleotide diveristy values below 0.0001 as '0'
zero_fold_pi_sect_primula$pi_pmega <- ifelse(zero_fold_pi_sect_primula$pi_pmega<=0.0005,0,zero_fold_pi_sect_primula$pi_pmega)
zero_fold_pi_sect_primula$pi_preni <- ifelse(zero_fold_pi_sect_primula$pi_preni<=0.0005,0,zero_fold_pi_sect_primula$pi_preni)
zero_fold_pi_sect_primula$pi_pveris <- ifelse(zero_fold_pi_sect_primula$pi_pveris<=0.0005,0,zero_fold_pi_sect_primula$pi_pveris)
zero_fold_pi_sect_primula$pi_pvulg <- ifelse(zero_fold_pi_sect_primula$pi_pvulg<=0.0005,0,zero_fold_pi_sect_primula$pi_pvulg)
zero_fold_pi_sect_primula$pi_pgrand_B <- ifelse(zero_fold_pi_sect_primula$pi_pgrand_B<=0.0005,0,zero_fold_pi_sect_primula$pi_pgrand_B)

four_fold_pi_grandis_A$pi_A <- ifelse(four_fold_pi_grandis_A$pi_A<=0.0005,0,four_fold_pi_grandis_A$pi_A)
zero_fold_pi_grandis_A$pi_A <- ifelse(zero_fold_pi_grandis_A$pi_A<=0.0005,0,zero_fold_pi_grandis_A$pi_A)
zero_fold_pi_grandis_A$pi_A <- ifelse(zero_fold_pi_grandis_A$pi_A<=0,NA,zero_fold_pi_grandis_A$pi_A)

four_fold_pi_grandis_B$pi_A <- ifelse(four_fold_pi_grandis_B$pi_A<=0.0005,0,four_fold_pi_grandis_B$pi_A)
zero_fold_pi_grandis_B$pi_A <- ifelse(zero_fold_pi_grandis_B$pi_A<=0.0005,0,zero_fold_pi_grandis_B$pi_A)
zero_fold_pi_grandis_B$pi_A <- ifelse(zero_fold_pi_grandis_B$pi_A<=0,NA,zero_fold_pi_grandis_B$pi_A)

piN_piS_pela <- zero_fold_pi_sect_primula %>% mutate(piN_piS_pela = zero_fold_pi_sect_primula$pi_pelatior/four_fold_pi_sect_primula$pi_pelatior)
piN_piS_pmega <- zero_fold_pi_sect_primula %>% mutate(piN_piS_pmega = zero_fold_pi_sect_primula$pi_pmega/four_fold_pi_sect_primula$pi_pmega)
piN_piS_preni <- zero_fold_pi_sect_primula %>% mutate(piN_piS_preni = zero_fold_pi_sect_primula$pi_preni/four_fold_pi_sect_primula$pi_preni)
piN_piS_pveris <- zero_fold_pi_sect_primula %>% mutate(piN_piS_pveris = zero_fold_pi_sect_primula$pi_pveris/four_fold_pi_sect_primula$pi_pveris)
piN_piS_pvulg <- zero_fold_pi_sect_primula %>% mutate(piN_piS_pvulg = zero_fold_pi_sect_primula$pi_pvulg/four_fold_pi_sect_primula$pi_pvulg)
piN_piS_pgrandis <- zero_fold_pi_sect_primula %>% mutate(piN_piS_grandis = zero_fold_pi_sect_primula$pi_pgrand/four_fold_pi_sect_primula$pi_pgrand)
piN_piS_subgenome_A <- zero_fold_pi_grandis_A %>% mutate(piN_piS_sub_A = zero_fold_pi_grandis_A$pi_A/four_fold_pi_grandis_A$pi_A)
piN_piS_subgenome_B <- zero_fold_pi_grandis_B %>% mutate(piN_piS_sub_B = zero_fold_pi_grandis_B$pi_A/four_fold_pi_grandis_B$pi_A)

is.na(piN_piS_pela) <- sapply(piN_piS_pela, is.infinite)
is.na(piN_piS_pmega) <- sapply(piN_piS_pmega, is.infinite)
is.na(piN_piS_preni) <- sapply(piN_piS_preni, is.infinite)
is.na(piN_piS_pveris) <- sapply(piN_piS_pveris, is.infinite)
is.na(piN_piS_pvulg) <- sapply(piN_piS_pvulg, is.infinite)
is.na(piN_piS_pgrandis) <- sapply(piN_piS_pgrandis, is.infinite)
is.na(piN_piS_subgenome_A) <- sapply(piN_piS_subgenome_A, is.infinite)
is.na(piN_piS_subgenome_B) <- sapply(piN_piS_subgenome_B, is.infinite)

summary(piN_piS_pela$piN_piS_pela)
summary(piN_piS_pmega$piN_piS_pmega)
summary(piN_piS_preni$piN_piS_preni)
summary(piN_piS_pveris$piN_piS_pveris)
summary(piN_piS_pvulg$piN_piS_pvulg)
summary(piN_piS_pgrandis$piN_piS_grandis)
summary(piN_piS_subgenome_A$piN_piS_sub_A)
summary(piN_piS_subgenome_B$piN_piS_sub_B)

############################
#PLOT piN / piS
############################

piN_piS_pela_melted <- melt(piN_piS_pela, id='windowID')  %>%  filter(variable %in% c('piN_piS_pela'))
piN_piS_pmega_melted <- melt(piN_piS_pmega, id='windowID')  %>%  filter(variable %in% c('piN_piS_pmega'))
piN_piS_preni_melted <- melt(piN_piS_preni, id='windowID')  %>%  filter(variable %in% c('piN_piS_preni'))
piN_piS_pveris_melted <- melt(piN_piS_pveris, id='windowID')  %>%  filter(variable %in% c('piN_piS_pveris'))
piN_piS_pvulg_melted <- melt(piN_piS_pvulg, id='windowID')  %>%  filter(variable %in% c('piN_piS_pvulg'))
#piN_piS_pgrand_melted <- melt(piN_piS_pgrand, id='windowID')  %>%  filter(variable %in% c('piN_piS_pgrand'))
piN_piS_subgenome_A_melted <- melt(piN_piS_subgenome_A, id='windowID')  %>%  filter(variable %in% c('piN_piS_sub_A'))
piN_piS_subgenome_B_melted <- melt(piN_piS_subgenome_B, id='windowID')  %>%  filter(variable %in% c('piN_piS_sub_B'))

df_list <- list(piN_piS_pela_melted, piN_piS_pmega_melted, piN_piS_preni_melted, piN_piS_pveris_melted, piN_piS_pvulg_melted,
                piN_piS_subgenome_A_melted,piN_piS_subgenome_B_melted)
piN_piS_sect_primula_melted <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)
colnames(piN_piS_sect_primula_melted) <- c('windowID','species','piN_piS')
piN_piS_sect_primula_melted$piN_piS <-as.numeric(piN_piS_sect_primula_melted$piN_piS)

#Plot as violin plot

ggplot(piN_piS_sect_primula_melted, aes(x=species, y=piN_piS, fill=species)) + 
  geom_violin(trim = T, na.rm = T, scale = "width", colour = "white") + 
  geom_boxplot(width=0.30, fill="white", na.rm = T, outlier.colour = "#565656", outlier.size = 0.6) + 
  ylim(0,5) + theme_classic(base_size = 11) +
  labs(title="", x="", y = "pi-Non-Synonymous / pi-Synonymous") + 
  scale_fill_manual(values=c("#ff7f00", "#75CCCB", "#FFE733","#4daf4a", "#f781bf", "#C73333","#377eb8"), 
                    labels=c("P. elatior","P. megaseifolia","P. renifolia","P. veris","P. vulgaris","P. grandis_A","P. grandis_B")) +
  scale_x_discrete(labels=c("P. elatior","P. megaseifolia","P. renifolia","P. veris","P. vulgaris","P. grandis_A","P. grandis_B")) +
  theme(axis.text.x = element_text(face = "italic", size=11, color = "black")) +
  theme(axis.text.y = element_text(size=10, color = "black"))

#Plot as ridgeline plot

ggplot(piN_piS_sect_primula_melted, aes(x = piN_piS, y = species, fill=species)) +
  geom_density_ridges(scale = 5, alpha = 0.75, bandwidth = 0.5, quantile_lines = TRUE, quantile_fun = median) + 
  scale_x_continuous(limits = c(0,2)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_manual(values=c("#ff7f00", "#75CCCB", "#FFE733","#4daf4a", "#f781bf", "#C73333","#377eb8"), 
                    labels=c("P. elatior","P. megaseifolia","P. renifolia","P. veris","P. vulgaris","P. grandis_A","P. grandis_B")) +
  scale_y_discrete(labels=c("P. elatior","P. megaseifolia","P. renifolia","P. veris","P. vulgaris","P. grandis_A","P. grandis_B")) +
  theme(axis.text.x = element_text(size=11, color = "black")) +
  theme(axis.text.y = element_text(face = "italic", size=10, color = "black")) +
  labs(title = "Distribution of pi-Non-Synonymous / pi-Synonymous", x = "pi-Non-Synonymous / pi-Synonymous") 

######################################################################################################
#Results from Simon Martin's scripts "popFreq"
######################################################################################################

four_fold_tajima_sect_primula <- read.csv("4fold.popgenWindows.Tajima.sect_primula.csv", header = T)
four_fold_tajima_grandis_A <- read.csv("4fold_pgrandis_A.popgenWindows.Tajima.csv", header = T)
four_fold_tajima_grandis_B <- read.csv("4fold_pgrandis_B.popgenWindows.Tajima.csv", header = T)

four_fold_tajima_sect_primula_subgen_A <- read.csv("4fold.popgenWindows.Tajima.sect_primula.a.csv", header = T)

##########################################################
#Segregating sites
sum(four_fold_tajima_sect_primula[, c("S_pelatior")], na.rm=TRUE)
sum(four_fold_tajima_sect_primula[, c("S_pmega")], na.rm=TRUE)
sum(four_fold_tajima_sect_primula[, c("S_preni")], na.rm=TRUE)
sum(four_fold_tajima_sect_primula[, c("S_pveris")], na.rm=TRUE)
sum(four_fold_tajima_sect_primula[, c("S_pvulg")], na.rm=TRUE)
#sum(four_fold_tajima_sect_primula[, c("S_pgrand_B")], na.rm=TRUE)
sum(four_fold_tajima_grandis_A[, c("S_A")], na.rm=TRUE)
sum(four_fold_tajima_grandis_B[, c("S_A")], na.rm=TRUE)

sum(four_fold_tajima_grandis_A[, c("sites")], na.rm=TRUE)
sum(four_fold_tajima_grandis_B[, c("sites")], na.rm=TRUE)

sum(four_fold_tajima_sect_primula_subgen_A[, c("S_pgrand_A")], na.rm=TRUE)

###########################################################
#Watterson's theta
summary(four_fold_tajima_sect_primula$thetaW_pelatior)
summary(four_fold_tajima_sect_primula$thetaW_pmega)
summary(four_fold_tajima_sect_primula$thetaW_preni)
summary(four_fold_tajima_sect_primula$thetaW_pveris)
summary(four_fold_tajima_sect_primula$thetaW_pvulg)
#summary(four_fold_tajima_sect_primula$thetaW_pgrand_B)
summary(four_fold_tajima_grandis_A$thetaW_A)
summary(four_fold_tajima_grandis_B$thetaW_A)

summary(four_fold_tajima_sect_primula_subgen_A$thetaW_pgrand_A)

#############################################################
#Theta Pi
summary(four_fold_tajima_sect_primula$thetaPi_pelatior)
summary(four_fold_tajima_sect_primula$thetaPi_pmega)
summary(four_fold_tajima_sect_primula$thetaPi_preni)
summary(four_fold_tajima_sect_primula$thetaPi_pveris)
summary(four_fold_tajima_sect_primula$thetaPi_pvulg)
summary(four_fold_tajima_sect_primula$thetaPi_pgrand_B)
summary(four_fold_tajima_grandis_A$thetaPi_A)
summary(four_fold_tajima_grandis_B$thetaPi_A)

summary(four_fold_tajima_sect_primula_subgen_A$thetaPi_pgrand_A)

###########################################################
#Tajima's D
###########################################################

summary(four_fold_tajima_sect_primula$TajD_pelatior)
summary(four_fold_tajima_sect_primula$TajD_pmega)
summary(four_fold_tajima_sect_primula$TajD_preni)
summary(four_fold_tajima_sect_primula$TajD_pveris)
summary(four_fold_tajima_sect_primula$TajD_pvulg)
#summary(four_fold_tajima_sect_primula$TajD_pgrand_B)
summary(four_fold_tajima_grandis_A$TajD_A)
summary(four_fold_tajima_grandis_B$TajD_A)

summary(four_fold_tajima_sect_primula_subgen_A$TajD_pgrand_A)

four_fold_tajima_sect_primula_melted <- melt(four_fold_tajima_sect_primula, id='windowID')  %>%  
  filter(variable %in% c('TajD_pelatior','TajD_pmega','TajD_preni','TajD_pveris','TajD_pvulg'))

colnames(four_fold_tajima_grandis_A)[11] <- "TajD_pgrandis_A" #change name of column 11
four_fold_tajima_grandis_A_melted <- melt(four_fold_tajima_grandis_A,id='windowID')  %>%  
  filter(variable == 'TajD_pgrandis_A')

colnames(four_fold_tajima_grandis_B)[11] <- "TajD_pgrandis_B" #change name of column 11
four_fold_tajima_grandis_B_melted <- melt(four_fold_tajima_grandis_B,id='windowID')  %>%  
  filter(variable == 'TajD_pgrandis_B')

#merge the three dataframes (sect_primula, grandis A and grandis B)

df_list <- list(four_fold_tajima_sect_primula_melted, four_fold_tajima_grandis_A_melted, four_fold_tajima_grandis_B_melted)
four_fold_tajima_sect_primula_melted_2 <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)

colnames(four_fold_tajima_sect_primula_melted_2) <- c('windowID','species','Taj_D')
four_fold_tajima_sect_primula_melted_2$Taj_D <-as.numeric(four_fold_tajima_sect_primula_melted_2$Taj_D)

ggplot(four_fold_tajima_sect_primula_melted_2, aes(x = Taj_D, y = species, fill=species)) +
  geom_density_ridges(scale = 5, alpha = 0.75, bandwidth = 0.25, quantile_lines = TRUE, quantile_fun = median) + 
  scale_x_continuous(limits = c(-3,3)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_manual(values=c("#ff7f00", "#75CCCB", "#FFE733","#4daf4a", "#f781bf", "#C73333","#377eb8"), 
                    labels=c("P. elatior","P. megaseifolia","P. renifolia","P. veris","P. vulgaris","P. grandis_A","P. grandis_B")) +
  scale_y_discrete(labels=c("P. elatior","P. megaseifolia","P. renifolia","P. veris","P. vulgaris","P. grandis_A","P. grandis_B")) +
  theme(axis.text.x = element_text(size=11, color = "black")) +
  theme(axis.text.y = element_text(face = "italic", size=10, color = "black")) +
  labs(title = "Distribution of Tajima's D", x = "Tajima's D")
```
