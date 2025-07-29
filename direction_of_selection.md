# Direction of Selection (DoS)
Finally, we estimated strength and Direction of Selection (DoS) of all genes across the genome based on the McDonald-Kreitman (MK) test using the script degenotate.py (https://github.com/harvardinformatics/degenotate). This approach allowed us to estimate the strength of positive and purifying selection of each gene and determine whether there are differences among the polyploid, selfing P. grandis and the diploid, heterostylous species of Primula sect. Primula. For this analysis, six individuals of P. veris were used as outgroup.
```sh
VCF=/home/ubuntu/emiliano/4th_chp/4_variant_calling/other_primula/subgenome_A/filtering/merge_grandis/merged_allpops_OGs.a.vcf.gz
REF=/home/ubuntu/emiliano/p_grandis_ref_genome/out_JBAT_review4.FINAL_rename.fasta
ANN=/home/ubuntu/giacomo-3/convEvo/mcscan_convEvo/grandis/pgra.allGenes.gtf

source /home/ubuntu/etienne/progs/miniconda3/bin/activate emiliano3.10
python /home/ubuntu/etienne/progs/miniconda3/envs/emiliano3.10/bin/degenotate.py \
-a $ANN \
-g $REF \
-v $VCF \
-e EXCLUDE_IND_SUB_A.txt \
-u OUTRGROUPS_A.txt \
-o sub_A \
--overwrite

VCF=/home/ubuntu/emiliano/4th_chp/4_variant_calling/other_primula/subgenome_B/filtering/merge_grandis/merged_allpops_OGs.b.vcf.gz
source /home/ubuntu/etienne/progs/miniconda3/bin/activate emiliano3.10
python /home/ubuntu/etienne/progs/miniconda3/envs/emiliano3.10/bin/degenotate.py \
-a $ANN \
-g $REF \
-v $VCF \
-e EXCLUDE_IND_SUB_B.txt \
-u OUTRGROUPS_B.txt \
-o sub_B \
--overwrite
```
Content of EXCLUDE_IND_SUB_A.txt
```
Pver_ver_WE02.SubGenA.bam
Pela_ela_WG05.SubGenA.bam
Pela_ela_WE09.SubGenA.bam
Pela_palla_WD02.SubGenA.bam
Pela_palla_WC03.SubGenA.bam
Pela_pseudo_WA11.SubGenA.bam
Pela_pseudo_WH11.SubGenA.bam
Pela_pseudo_WG12.SubGenA.bam
Pela_intri_WG11.SubGenA.bam
Pela_leuco_WD01.SubGenA.bam
Pela_leuco_WC02.SubGenA.bam
Pela_leuco_WC03.SubGenA.bam
Pela_leuco_WB04.SubGenA.bam
Pela_mey_WE08.SubGenA.bam
Pela_mey_WE09.SubGenA.bam
Pela_mey_WE10.SubGenA.bam
Pela_mey_WF11.SubGenA.bam
Pela_cordi_WD08.SubGenA.bam
Pela_cordi_WD09.SubGenA.bam
Pela_cordi_WD10.SubGenA.bam
Pela_cordu_WE11.SubGenA.bam
Pela_palla_WH04.SubGenA.bam
Pela_palla_WA04.SubGenA.bam
Pela_palla_WB02.SubGenA.bam
Pela_palla_WA03.SubGenA.bam
Pela_palla_WB01.SubGenA.bam
Pela_mey_WF05.SubGenA.bam
Pela_mey_WD07.SubGenA.bam
Pjul_WH10.SubGenA.bam
Pjul_WG11.SubGenA.bam
Pjul_WF12.SubGenA.bam
Pjul_WE01.SubGenA.bam
Pmega_WC06.SubGenA.bam
Pmega_WC07.SubGenA.bam
Pmega_WF02.SubGenA.bam
Pmega_WF09.SubGenA.bam
Pmega_WF10.SubGenA.bam
Pmega_WE11.SubGenA.bam
Pmega_WD12.SubGenA.bam
Pmega_WD01.SubGenA.bam
Pmega_WD02.SubGenA.bam
Pmega_WD03.SubGenA.bam
Preni_WH01.SubGenA.bam
Preni_WH02.SubGenA.bam
Preni_WH05.SubGenA.bam
Preni_WA05.SubGenA.bam
Preni_WA06.SubGenA.bam
Preni_WA07.SubGenA.bam
Preni_WA08.SubGenA.bam
Preni_WA09.SubGenA.bam
Preni_WA10.SubGenA.bam
Preni_WA11.SubGenA.bam
Pver_col_WE06.SubGenA.bam
Pver_col_WH03.SubGenA.bam
Pver_col_WD09.SubGenA.bam
Pver_macro_WF06.SubGenA.bam
Pver_macro_WE07.SubGenA.bam
Pver_ver_WE12.SubGenA.bam
Pvul_vul_WF05.SubGenA.bam
Pvul_hetero_WA12.SubGenA.bam
Pvul_hetero_WA01.SubGenA.bam
Pvul_hetero_WA02.SubGenA.bam
Pvul_hetero_WA03.SubGenA.bam
Pvul_vul_WG02.SubGenA.bam
Pvul_vul_WG09.SubGenA.bam
Pvul_vul_WG10.SubGenA.bam
Pvul_vul_WF11.SubGenA.bam
Pvul_vul_WE12.SubGenA.bam
Pvul_vul_WE01.SubGenA.bam
Pvul_vul_WE02.SubGenA.bam
Pvul_vul_WE03.SubGenA.bam
Pvul_vul_WE04.SubGenA.bam
Pvul_vul_WE05.SubGenA.bam
Pvul_vul_WE06.SubGenA.bam
Pvul_sib_WC02.SubGenA.bam
Pvul_sib_WC03.SubGenA.bam
Pvul_sib_WC04.SubGenA.bam
Pvul_sib_WC05.SubGenA.bam
Pvul_sib_WC01.SubGenA.bam
Pvul_sib_WB03.SubGenA.bam
Psiamensis_WA09.SubGenA.bam
Psouliei_WB09.SubGenA.bam
Pbouveana_WC09.SubGenA.bam
Pbouveana_WD09.SubGenA.bam
Dionysia_WE08.SubGenA.bam
Pmistassinica_WE09.SubGenA.bam
Pverticilliata_WF08.SubGenA.bam
Pnutans_WF09.SubGenA.bam
```
EXCLUDE_IND_SUB_B.txt
```
Pver_ver_WE02.SubGenB.bam
Pela_ela_WG05.SubGenB.bam
Pela_ela_WE09.SubGenB.bam
Pela_palla_WD02.SubGenB.bam
Pela_palla_WC03.SubGenB.bam
Pela_pseudo_WA11.SubGenB.bam
Pela_pseudo_WH11.SubGenB.bam
Pela_pseudo_WG12.SubGenB.bam
Pela_intri_WG11.SubGenB.bam
Pela_leuco_WD01.SubGenB.bam
Pela_leuco_WC02.SubGenB.bam
Pela_leuco_WC03.SubGenB.bam
Pela_leuco_WB04.SubGenB.bam
Pela_mey_WE08.SubGenB.bam
Pela_mey_WE09.SubGenB.bam
Pela_mey_WE10.SubGenB.bam
Pela_mey_WF11.SubGenB.bam
Pela_cordi_WD08.SubGenB.bam
Pela_cordi_WD09.SubGenB.bam
Pela_cordi_WD10.SubGenB.bam
Pela_cordu_WE11.SubGenB.bam
Pela_palla_WH04.SubGenB.bam
Pela_palla_WA04.SubGenB.bam
Pela_palla_WB02.SubGenB.bam
Pela_palla_WA03.SubGenB.bam
Pela_palla_WB01.SubGenB.bam
Pela_mey_WF05.SubGenB.bam
Pela_mey_WD07.SubGenB.bam
Pjul_WH10.SubGenB.bam
Pjul_WG11.SubGenB.bam
Pjul_WF12.SubGenB.bam
Pjul_WE01.SubGenB.bam
Pmega_WC06.SubGenB.bam
Pmega_WC07.SubGenB.bam
Pmega_WF02.SubGenB.bam
Pmega_WF09.SubGenB.bam
Pmega_WF10.SubGenB.bam
Pmega_WE11.SubGenB.bam
Pmega_WD12.SubGenB.bam
Pmega_WD01.SubGenB.bam
Pmega_WD02.SubGenB.bam
Pmega_WD03.SubGenB.bam
Preni_WH01.SubGenB.bam
Preni_WH02.SubGenB.bam
Preni_WH05.SubGenB.bam
Preni_WA05.SubGenB.bam
Preni_WA06.SubGenB.bam
Preni_WA07.SubGenB.bam
Preni_WA08.SubGenB.bam
Preni_WA09.SubGenB.bam
Preni_WA10.SubGenB.bam
Preni_WA11.SubGenB.bam
Pver_col_WE06.SubGenB.bam
Pver_col_WH03.SubGenB.bam
Pver_col_WD09.SubGenB.bam
Pver_macro_WF06.SubGenB.bam
Pver_macro_WE07.SubGenB.bam
Pver_ver_WE12.SubGenB.bam
Pvul_vul_WF05.SubGenB.bam
Pvul_hetero_WA12.SubGenB.bam
Pvul_hetero_WA01.SubGenB.bam
Pvul_hetero_WA02.SubGenB.bam
Pvul_hetero_WA03.SubGenB.bam
Pvul_vul_WG02.SubGenB.bam
Pvul_vul_WG09.SubGenB.bam
Pvul_vul_WG10.SubGenB.bam
Pvul_vul_WF11.SubGenB.bam
Pvul_vul_WE12.SubGenB.bam
Pvul_vul_WE01.SubGenB.bam
Pvul_vul_WE02.SubGenB.bam
Pvul_vul_WE03.SubGenB.bam
Pvul_vul_WE04.SubGenB.bam
Pvul_vul_WE05.SubGenB.bam
Pvul_vul_WE06.SubGenB.bam
Pvul_sib_WC02.SubGenB.bam
Pvul_sib_WC03.SubGenB.bam
Pvul_sib_WC04.SubGenB.bam
Pvul_sib_WC05.SubGenB.bam
Pvul_sib_WC01.SubGenB.bam
Pvul_sib_WB03.SubGenB.bam
Psiamensis_WA09.SubGenB.bam
Psouliei_WB09.SubGenB.bam
Pbouveana_WC09.SubGenB.bam
Pbouveana_WD09.SubGenB.bam
Dionysia_WE08.SubGenB.bam
Pmistassinica_WE09.SubGenB.bam
Pverticilliata_WF08.SubGenB.bam
Pnutans_WF09.SubGenB.bam
```
Content of OUTRGROUPS_A.txt
```
Pver_macro_WA05.SubGenA.bam
Pver_macro_WH05.SubGenA.bam
Pver_macro_WG06.SubGenA.bam
Pver_macro_WF07.SubGenA.bam
Pver_macro_WH04.SubGenA.bam
Pver_macro_WG05.SubGenA.bam
```
Content of OUTRGROUPS_B.txt
```
Pver_macro_WA05.SubGenB.bam
Pver_macro_WH05.SubGenB.bam
Pver_macro_WG06.SubGenB.bam
Pver_macro_WF07.SubGenB.bam
Pver_macro_WH04.SubGenB.bam
Pver_macro_WG05.SubGenB.bam
```
R volcano plots
```R
library(tidyverse)
library(ggpubr)
library(ggplot2)

setwd("~/r_analysis/primula_grandis/databases/DoS")

#Code based on: https://github.com/sjswuitchik/compPopGen_ms/blob/master/Mirchandani_2023/mk/mk.Rmd

###################################################
#SUBGENOME A vs SUBGENOME B (using P. veris as an outgroup)
###################################################

mk_subA_OG_veris <- read_tsv("mk_sub_A_veris_5.tsv")
mk_subB_OG_veris <- read_tsv("mk_sub_B_veris_5.tsv")

tsinfo <- read_tsv("transcript-counts.tsv")

merge_subA <- full_join(mk_subA_OG_veris, tsinfo, by="transcript") %>% filter(is_longest=="yes") %>% select(gene, cds_length, pN, pS, dN, dS, f0, f2, f3, f4, pval, dos)
merge_subA_complete <- merge_subA

merge_subB <- full_join(mk_subB_OG_veris, tsinfo, by="transcript") %>% filter(is_longest=="yes") %>% select(gene, cds_length, pN, pS, dN, dS, f0, f2, f3, f4, pval, dos)
merge_subB_complete <- merge_subB

#Filters

merge_subAfiltered <- merge_subA_complete %>% mutate(logP = -1 * log10(pval), direction = ifelse(dos < 0, "neg", "pos"), sig = ifelse(logP>= 2, "yes", "no"), 
                              result = paste0(direction, ":", sig)) %>% # ADDS ONE COLUMN WITH THE logP VALUE AND THREE CONDITIONALS COLUMNS BASED ON THE SIGNIFICANCE 
  filter(!is.na(dos), cds_length <= 10000, !is.na(logP), dN+dS > 2, pN+pS > 2) %>%  #FILTERS SHORT TRANSCRIPTS, NA's and some dN,dS AND pNpS threshold
  mutate(gene_label = str_remove(gene, "gene-")) %>% #CREATES GENE LABEL
  mutate(gene_label = ifelse(str_detect(gene_label, "LOC"), "", gene_label))

merge_subB_filtered <- merge_subB_complete %>% mutate(logP = -1 * log10(pval), direction = ifelse(dos < 0, "neg", "pos"), sig = ifelse(logP>= 2, "yes", "no"), 
                                                        result = paste0(direction, ":", sig)) %>% # ADDS ONE COLUMN WITH THE logP VALUE AND THREE CONDITIONALS COLUMNS BASED ON THE SIGNIFICANCE 
  filter(!is.na(dos), cds_length <= 10000, !is.na(logP), dN+dS > 2, pN+pS > 2) %>%  #FILTERS SHORT TRANSCRIPTS, NA's and some dN,dS AND pNpS threshold
  mutate(gene_label = str_remove(gene, "gene-")) %>% #CREATES GENE LABEL
  mutate(gene_label = ifelse(str_detect(gene_label, "LOC"), "", gene_label))

#Plot

plot_volcano_A_veris <-ggplot(data=merge_subAfiltered, aes(y=logP, x=dos, color=result, alpha=sig)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values=c("grey", "black",  "grey", "black"), guide="none") +
  scale_alpha_discrete(range = c(0.6, 1), guide="none") +
  xlab("Direction of Selection (DoS)") + scale_x_continuous(limits = c(-1,1)) +
  ylab("-log10 P-value (MK test)") + scale_y_continuous(limits = c(0,3.5)) +
  geom_text(aes(label=ifelse(logP>2,as.character(gene_label),'')), hjust=-0.1,vjust=0.4,color="black") #ADDS LABEL IF SOME TRANSCRIPTS HAVE VERY HIGH PURYFING SELECTION

plot_volcano_B_veris <-ggplot(data=merge_subB_filtered, aes(y=logP, x=dos, color=result, alpha=sig)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values=c("grey", "black",  "grey", "black"), guide="none") +
  scale_alpha_discrete(range = c(0.6, 1), guide="none") +
#  scale_color_manual(values=c("#1f78b4", "#1f78b4",  "#e31a1c", "#e31a1c"), guide="none") +
#  scale_alpha_discrete(range = c(0.2, 0.8), guide="none") +
  xlab("Direction of Selection (DoS)") + scale_x_continuous(limits = c(-1,1)) +
  ylab("-log10 P-value (MK test)") + scale_y_continuous(limits = c(0,3.5)) +
  geom_text(aes(label=ifelse(logP>2,as.character(gene_label),'')), hjust=-0.1,vjust=0.4,color="black") 

ggarrange(plot_volcano_A_veris, plot_volcano_B_veris, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

##
```sh
VCF=/home/ubuntu/emiliano/4th_chp/4_variant_calling/other_primula/subgenome_B/filtering/merge_grandis/merged_allpops_OGs.b.vcf.gz
REF=/home/ubuntu/emiliano/p_grandis_ref_genome/out_JBAT_review4.FINAL_rename.fasta
ANN=/home/ubuntu/giacomo-3/convEvo/mcscan_convEvo/grandis/pgra.allGenes.gtf

SPECIES=(Pela_palla Pmega Preni Pvul_vul Pver_macro)

for i in ${SPECIES[@]}; do
source /home/ubuntu/etienne/progs/miniconda3/bin/activate emiliano3.10
python /home/ubuntu/etienne/progs/miniconda3/envs/emiliano3.10/bin/degenotate.py \
-a $ANN \
-g $REF \
-v $VCF \
-e EXCLUDE_FILES/EXCLUDE.$i.txt \
-u EXCLUDE_FILES/OUTRGROUPS_B.txt \
-o $i \
--overwrite

mv $i/mk.tsv $i/mk_$i.tsv; done
```
Content of all files in the script above:
```
::::::::::::::
EXCLUDE_FILES/EXCLUDE.Pela_palla.txt
::::::::::::::
Pela_ela_WG05.SubGenB.bam
Pela_ela_WE09.SubGenB.bam
Pela_pseudo_WA11.SubGenB.bam
Pela_pseudo_WH11.SubGenB.bam
Pela_pseudo_WG12.SubGenB.bam
Pela_intri_WG11.SubGenB.bam
Pela_leuco_WD01.SubGenB.bam
Pela_leuco_WC02.SubGenB.bam
Pela_leuco_WC03.SubGenB.bam
Pela_leuco_WB04.SubGenB.bam
Pela_mey_WE08.SubGenB.bam
Pela_mey_WE09.SubGenB.bam
Pela_mey_WE10.SubGenB.bam
Pela_mey_WF11.SubGenB.bam
Pela_cordi_WD08.SubGenB.bam
Pela_cordi_WD09.SubGenB.bam
Pela_cordi_WD10.SubGenB.bam
Pela_cordu_WE11.SubGenB.bam
Pela_palla_WA04.SubGenB.bam
Pela_mey_WF05.SubGenB.bam
Pela_mey_WD07.SubGenB.bam
Pjul_WH10.SubGenB.bam
Pjul_WG11.SubGenB.bam
Pjul_WF12.SubGenB.bam
Pjul_WE01.SubGenB.bam
Pmega_WC06.SubGenB.bam
Pmega_WC07.SubGenB.bam
Pmega_WF02.SubGenB.bam
Pmega_WF09.SubGenB.bam
Pmega_WF10.SubGenB.bam
Pmega_WE11.SubGenB.bam
Pmega_WD12.SubGenB.bam
Pmega_WD01.SubGenB.bam
Pmega_WD02.SubGenB.bam
Pmega_WD03.SubGenB.bam
Preni_WH01.SubGenB.bam
Preni_WH02.SubGenB.bam
Preni_WH05.SubGenB.bam
Preni_WA05.SubGenB.bam
Preni_WA06.SubGenB.bam
Preni_WA07.SubGenB.bam
Preni_WA08.SubGenB.bam
Preni_WA09.SubGenB.bam
Preni_WA10.SubGenB.bam
Preni_WA11.SubGenB.bam
Pver_col_WE06.SubGenB.bam
Pver_ver_WE02.SubGenB.bam
Pver_col_WH03.SubGenB.bam
Pver_col_WD09.SubGenB.bam
Pver_macro_WF06.SubGenB.bam
Pver_macro_WE07.SubGenB.bam
Pver_ver_WE12.SubGenB.bam
Pvul_vul_WF05.SubGenB.bam
Pvul_hetero_WA12.SubGenB.bam
Pvul_hetero_WA01.SubGenB.bam
Pvul_hetero_WA02.SubGenB.bam
Pvul_hetero_WA03.SubGenB.bam
Pvul_vul_WG02.SubGenB.bam
Pvul_vul_WG09.SubGenB.bam
Pvul_vul_WG10.SubGenB.bam
Pvul_vul_WF11.SubGenB.bam
Pvul_vul_WE12.SubGenB.bam
Pvul_vul_WE01.SubGenB.bam
Pvul_vul_WE02.SubGenB.bam
Pvul_vul_WE03.SubGenB.bam
Pvul_vul_WE04.SubGenB.bam
Pvul_vul_WE05.SubGenB.bam
Pvul_vul_WE06.SubGenB.bam
Pvul_sib_WC02.SubGenB.bam
Pvul_sib_WC03.SubGenB.bam
Pvul_sib_WC04.SubGenB.bam
Pvul_sib_WC05.SubGenB.bam
Pvul_sib_WC01.SubGenB.bam
Pvul_sib_WB03.SubGenB.bam
grandis_WA05.bam
grandis_WB03.bam
grandis_WB04.bam
grandis_WB10.bam
grandis_WC02.bam
grandis_WC09.bam
grandis_WD08.bam
grandis_WE08.bam
grandis_WF03.bam
grandis_WG03.bam
::::::::::::::
EXCLUDE_FILES/EXCLUDE.Pmega.txt
::::::::::::::
Pela_ela_WG05.SubGenB.bam
Pela_ela_WE09.SubGenB.bam
Pela_palla_WD02.SubGenB.bam
Pela_palla_WC03.SubGenB.bam
Pela_pseudo_WA11.SubGenB.bam
Pela_pseudo_WH11.SubGenB.bam
Pela_pseudo_WG12.SubGenB.bam
Pela_intri_WG11.SubGenB.bam
Pela_leuco_WD01.SubGenB.bam
Pela_leuco_WC02.SubGenB.bam
Pela_leuco_WC03.SubGenB.bam
Pela_leuco_WB04.SubGenB.bam
Pela_mey_WE08.SubGenB.bam
Pela_mey_WE09.SubGenB.bam
Pela_mey_WE10.SubGenB.bam
Pela_mey_WF11.SubGenB.bam
Pela_cordi_WD08.SubGenB.bam
Pela_cordi_WD09.SubGenB.bam
Pela_cordi_WD10.SubGenB.bam
Pela_cordu_WE11.SubGenB.bam
Pela_palla_WH04.SubGenB.bam
Pela_palla_WA04.SubGenB.bam
Pela_palla_WB02.SubGenB.bam
Pela_palla_WA03.SubGenB.bam
Pela_palla_WB01.SubGenB.bam
Pela_mey_WF05.SubGenB.bam
Pela_mey_WD07.SubGenB.bam
Pjul_WH10.SubGenB.bam
Pjul_WG11.SubGenB.bam
Pjul_WF12.SubGenB.bam
Pjul_WE01.SubGenB.bam
Pmega_WF02.SubGenB.bam
Pmega_WF09.SubGenB.bam
Pmega_WF10.SubGenB.bam
Pmega_WD01.SubGenB.bam
Preni_WH01.SubGenB.bam
Preni_WH02.SubGenB.bam
Preni_WH05.SubGenB.bam
Preni_WA05.SubGenB.bam
Preni_WA06.SubGenB.bam
Preni_WA07.SubGenB.bam
Preni_WA08.SubGenB.bam
Preni_WA09.SubGenB.bam
Preni_WA10.SubGenB.bam
Preni_WA11.SubGenB.bam
Pver_col_WE06.SubGenB.bam
Pver_ver_WE02.SubGenB.bam
Pver_col_WH03.SubGenB.bam
Pver_col_WD09.SubGenB.bam
Pver_macro_WF06.SubGenB.bam
Pver_macro_WE07.SubGenB.bam
Pver_ver_WE12.SubGenB.bam
Pvul_vul_WF05.SubGenB.bam
Pvul_hetero_WA12.SubGenB.bam
Pvul_hetero_WA01.SubGenB.bam
Pvul_hetero_WA02.SubGenB.bam
Pvul_hetero_WA03.SubGenB.bam
Pvul_vul_WG02.SubGenB.bam
Pvul_vul_WG09.SubGenB.bam
Pvul_vul_WG10.SubGenB.bam
Pvul_vul_WF11.SubGenB.bam
Pvul_vul_WE12.SubGenB.bam
Pvul_vul_WE01.SubGenB.bam
Pvul_vul_WE02.SubGenB.bam
Pvul_vul_WE03.SubGenB.bam
Pvul_vul_WE04.SubGenB.bam
Pvul_vul_WE05.SubGenB.bam
Pvul_vul_WE06.SubGenB.bam
Pvul_sib_WC02.SubGenB.bam
Pvul_sib_WC03.SubGenB.bam
Pvul_sib_WC04.SubGenB.bam
Pvul_sib_WC05.SubGenB.bam
Pvul_sib_WC01.SubGenB.bam
Pvul_sib_WB03.SubGenB.bam
grandis_WA05.bam
grandis_WB03.bam
grandis_WB04.bam
grandis_WB10.bam
grandis_WC02.bam
grandis_WC09.bam
grandis_WD08.bam
grandis_WE08.bam
grandis_WF03.bam
grandis_WG03.bam
::::::::::::::
EXCLUDE_FILES/EXCLUDE.Preni.txt
::::::::::::::
Pela_ela_WG05.SubGenB.bam
Pela_ela_WE09.SubGenB.bam
Pela_palla_WD02.SubGenB.bam
Pela_palla_WC03.SubGenB.bam
Pela_pseudo_WA11.SubGenB.bam
Pela_pseudo_WH11.SubGenB.bam
Pela_pseudo_WG12.SubGenB.bam
Pela_intri_WG11.SubGenB.bam
Pela_leuco_WD01.SubGenB.bam
Pela_leuco_WC02.SubGenB.bam
Pela_leuco_WC03.SubGenB.bam
Pela_leuco_WB04.SubGenB.bam
Pela_mey_WE08.SubGenB.bam
Pela_mey_WE09.SubGenB.bam
Pela_mey_WE10.SubGenB.bam
Pela_mey_WF11.SubGenB.bam
Pela_cordi_WD08.SubGenB.bam
Pela_cordi_WD09.SubGenB.bam
Pela_cordi_WD10.SubGenB.bam
Pela_cordu_WE11.SubGenB.bam
Pela_palla_WH04.SubGenB.bam
Pela_palla_WA04.SubGenB.bam
Pela_palla_WB02.SubGenB.bam
Pela_palla_WA03.SubGenB.bam
Pela_palla_WB01.SubGenB.bam
Pela_mey_WF05.SubGenB.bam
Pela_mey_WD07.SubGenB.bam
Pjul_WH10.SubGenB.bam
Pjul_WG11.SubGenB.bam
Pjul_WF12.SubGenB.bam
Pjul_WE01.SubGenB.bam
Pmega_WC06.SubGenB.bam
Pmega_WC07.SubGenB.bam
Pmega_WF02.SubGenB.bam
Pmega_WF09.SubGenB.bam
Pmega_WF10.SubGenB.bam
Pmega_WE11.SubGenB.bam
Pmega_WD12.SubGenB.bam
Pmega_WD01.SubGenB.bam
Pmega_WD02.SubGenB.bam
Pmega_WD03.SubGenB.bam
Pver_col_WE06.SubGenB.bam
Pver_ver_WE02.SubGenB.bam
Pver_col_WH03.SubGenB.bam
Pver_col_WD09.SubGenB.bam
Pver_macro_WF06.SubGenB.bam
Pver_macro_WE07.SubGenB.bam
Pver_ver_WE12.SubGenB.bam
Pvul_vul_WF05.SubGenB.bam
Pvul_hetero_WA12.SubGenB.bam
Pvul_hetero_WA01.SubGenB.bam
Pvul_hetero_WA02.SubGenB.bam
Pvul_hetero_WA03.SubGenB.bam
Pvul_vul_WG02.SubGenB.bam
Pvul_vul_WG09.SubGenB.bam
Pvul_vul_WG10.SubGenB.bam
Pvul_vul_WF11.SubGenB.bam
Pvul_vul_WE12.SubGenB.bam
Pvul_vul_WE01.SubGenB.bam
Pvul_vul_WE02.SubGenB.bam
Pvul_vul_WE03.SubGenB.bam
Pvul_vul_WE04.SubGenB.bam
Pvul_vul_WE05.SubGenB.bam
Pvul_vul_WE06.SubGenB.bam
Pvul_sib_WC02.SubGenB.bam
Pvul_sib_WC03.SubGenB.bam
Pvul_sib_WC04.SubGenB.bam
Pvul_sib_WC05.SubGenB.bam
Pvul_sib_WC01.SubGenB.bam
Pvul_sib_WB03.SubGenB.bam
grandis_WA05.bam
grandis_WB03.bam
grandis_WB04.bam
grandis_WB10.bam
grandis_WC02.bam
grandis_WC09.bam
grandis_WD08.bam
grandis_WE08.bam
grandis_WF03.bam
grandis_WG03.bam
::::::::::::::
EXCLUDE_FILES/EXCLUDE.Pver_macro.txt
::::::::::::::
Pela_ela_WG05.SubGenB.bam
Pela_ela_WE09.SubGenB.bam
Pela_palla_WD02.SubGenB.bam
Pela_palla_WC03.SubGenB.bam
Pela_pseudo_WA11.SubGenB.bam
Pela_pseudo_WH11.SubGenB.bam
Pela_pseudo_WG12.SubGenB.bam
Pela_intri_WG11.SubGenB.bam
Pela_leuco_WD01.SubGenB.bam
Pela_leuco_WC02.SubGenB.bam
Pela_leuco_WC03.SubGenB.bam
Pela_leuco_WB04.SubGenB.bam
Pela_mey_WE08.SubGenB.bam
Pela_mey_WE09.SubGenB.bam
Pela_mey_WE10.SubGenB.bam
Pela_mey_WF11.SubGenB.bam
Pela_cordi_WD08.SubGenB.bam
Pela_cordi_WD09.SubGenB.bam
Pela_cordi_WD10.SubGenB.bam
Pela_cordu_WE11.SubGenB.bam
Pela_palla_WH04.SubGenB.bam
Pela_palla_WA04.SubGenB.bam
Pela_palla_WB02.SubGenB.bam
Pela_palla_WA03.SubGenB.bam
Pela_palla_WB01.SubGenB.bam
Pela_mey_WF05.SubGenB.bam
Pela_mey_WD07.SubGenB.bam
Pjul_WH10.SubGenB.bam
Pjul_WG11.SubGenB.bam
Pjul_WF12.SubGenB.bam
Pjul_WE01.SubGenB.bam
Pmega_WC06.SubGenB.bam
Pmega_WC07.SubGenB.bam
Pmega_WF02.SubGenB.bam
Pmega_WF09.SubGenB.bam
Pmega_WF10.SubGenB.bam
Pmega_WE11.SubGenB.bam
Pmega_WD12.SubGenB.bam
Pmega_WD01.SubGenB.bam
Pmega_WD02.SubGenB.bam
Pmega_WD03.SubGenB.bam
Preni_WH01.SubGenB.bam
Preni_WH02.SubGenB.bam
Preni_WH05.SubGenB.bam
Preni_WA05.SubGenB.bam
Preni_WA06.SubGenB.bam
Preni_WA07.SubGenB.bam
Preni_WA08.SubGenB.bam
Preni_WA09.SubGenB.bam
Preni_WA10.SubGenB.bam
Preni_WA11.SubGenB.bam
Pver_col_WE06.SubGenB.bam
Pver_ver_WE02.SubGenB.bam
Pver_col_WH03.SubGenB.bam
Pver_col_WD09.SubGenB.bam
Pver_ver_WE12.SubGenB.bam
Pvul_vul_WF05.SubGenB.bam
Pvul_hetero_WA12.SubGenB.bam
Pvul_hetero_WA01.SubGenB.bam
Pvul_hetero_WA02.SubGenB.bam
Pvul_hetero_WA03.SubGenB.bam
Pvul_vul_WE02.SubGenB.bam
Pvul_vul_WE03.SubGenB.bam
Pvul_vul_WE04.SubGenB.bam
Pvul_vul_WE05.SubGenB.bam
Pvul_vul_WE06.SubGenB.bam
Pvul_sib_WC02.SubGenB.bam
Pvul_sib_WC03.SubGenB.bam
Pvul_sib_WC04.SubGenB.bam
Pvul_sib_WC05.SubGenB.bam
Pvul_sib_WC01.SubGenB.bam
Pvul_sib_WB03.SubGenB.bam
grandis_WA05.bam
grandis_WB03.bam
grandis_WB04.bam
grandis_WB10.bam
grandis_WC02.bam
grandis_WC09.bam
grandis_WD08.bam
grandis_WE08.bam
grandis_WF03.bam
grandis_WG03.bam
::::::::::::::
EXCLUDE_FILES/EXCLUDE.Pvul_vul.txt
::::::::::::::
Pela_ela_WG05.SubGenB.bam
Pela_ela_WE09.SubGenB.bam
Pela_palla_WD02.SubGenB.bam
Pela_palla_WC03.SubGenB.bam
Pela_pseudo_WA11.SubGenB.bam
Pela_pseudo_WH11.SubGenB.bam
Pela_pseudo_WG12.SubGenB.bam
Pela_intri_WG11.SubGenB.bam
Pela_leuco_WD01.SubGenB.bam
Pela_leuco_WC02.SubGenB.bam
Pela_leuco_WC03.SubGenB.bam
Pela_leuco_WB04.SubGenB.bam
Pela_mey_WE08.SubGenB.bam
Pela_mey_WE09.SubGenB.bam
Pela_mey_WE10.SubGenB.bam
Pela_mey_WF11.SubGenB.bam
Pela_cordi_WD08.SubGenB.bam
Pela_cordi_WD09.SubGenB.bam
Pela_cordi_WD10.SubGenB.bam
Pela_cordu_WE11.SubGenB.bam
Pela_palla_WH04.SubGenB.bam
Pela_palla_WA04.SubGenB.bam
Pela_palla_WB02.SubGenB.bam
Pela_palla_WA03.SubGenB.bam
Pela_palla_WB01.SubGenB.bam
Pela_mey_WF05.SubGenB.bam
Pela_mey_WD07.SubGenB.bam
Pjul_WH10.SubGenB.bam
Pjul_WG11.SubGenB.bam
Pjul_WF12.SubGenB.bam
Pjul_WE01.SubGenB.bam
Pmega_WC06.SubGenB.bam
Pmega_WC07.SubGenB.bam
Pmega_WF02.SubGenB.bam
Pmega_WF09.SubGenB.bam
Pmega_WF10.SubGenB.bam
Pmega_WE11.SubGenB.bam
Pmega_WD12.SubGenB.bam
Pmega_WD01.SubGenB.bam
Pmega_WD02.SubGenB.bam
Pmega_WD03.SubGenB.bam
Preni_WH01.SubGenB.bam
Preni_WH02.SubGenB.bam
Preni_WH05.SubGenB.bam
Preni_WA05.SubGenB.bam
Preni_WA06.SubGenB.bam
Preni_WA07.SubGenB.bam
Preni_WA08.SubGenB.bam
Preni_WA09.SubGenB.bam
Preni_WA10.SubGenB.bam
Preni_WA11.SubGenB.bam
Pver_col_WE06.SubGenB.bam
Pver_ver_WE02.SubGenB.bam
Pver_col_WH03.SubGenB.bam
Pver_col_WD09.SubGenB.bam
Pver_macro_WF06.SubGenB.bam
Pver_macro_WE07.SubGenB.bam
Pver_ver_WE12.SubGenB.bam
Pvul_vul_WF05.SubGenB.bam
Pvul_hetero_WA12.SubGenB.bam
Pvul_hetero_WA01.SubGenB.bam
Pvul_hetero_WA02.SubGenB.bam
Pvul_hetero_WA03.SubGenB.bam
Pvul_sib_WC02.SubGenB.bam
Pvul_sib_WC03.SubGenB.bam
Pvul_sib_WC04.SubGenB.bam
Pvul_sib_WC05.SubGenB.bam
Pvul_sib_WC01.SubGenB.bam
Pvul_sib_WB03.SubGenB.bam
grandis_WA05.bam
grandis_WB03.bam
grandis_WB04.bam
grandis_WB10.bam
grandis_WC02.bam
grandis_WC09.bam
grandis_WD08.bam
grandis_WE08.bam
grandis_WF03.bam
grandis_WG03.bam
::::::::::::::
EXCLUDE_FILES/EXCLUDE.txt
::::::::::::::
Pela_ela_WG05.SubGenB.bam
Pela_ela_WE09.SubGenB.bam
Pela_palla_WD02.SubGenB.bam
Pela_palla_WC03.SubGenB.bam
Pela_pseudo_WA11.SubGenB.bam
Pela_pseudo_WH11.SubGenB.bam
Pela_pseudo_WG12.SubGenB.bam
Pela_intri_WG11.SubGenB.bam
Pela_leuco_WD01.SubGenB.bam
Pela_leuco_WC02.SubGenB.bam
Pela_leuco_WC03.SubGenB.bam
Pela_leuco_WB04.SubGenB.bam
Pela_mey_WE08.SubGenB.bam
Pela_mey_WE09.SubGenB.bam
Pela_mey_WE10.SubGenB.bam
Pela_mey_WF11.SubGenB.bam
Pela_cordi_WD08.SubGenB.bam
Pela_cordi_WD09.SubGenB.bam
Pela_cordi_WD10.SubGenB.bam
Pela_cordu_WE11.SubGenB.bam
Pela_palla_WH04.SubGenB.bam
Pela_palla_WA04.SubGenB.bam
Pela_palla_WB02.SubGenB.bam
Pela_palla_WA03.SubGenB.bam
Pela_palla_WB01.SubGenB.bam
Pela_mey_WF05.SubGenB.bam
Pela_mey_WD07.SubGenB.bam
Pjul_WH10.SubGenB.bam
Pjul_WG11.SubGenB.bam
Pjul_WF12.SubGenB.bam
Pjul_WE01.SubGenB.bam
Pmega_WC06.SubGenB.bam
Pmega_WC07.SubGenB.bam
Pmega_WF02.SubGenB.bam
Pmega_WF09.SubGenB.bam
Pmega_WF10.SubGenB.bam
Pmega_WE11.SubGenB.bam
Pmega_WD12.SubGenB.bam
Pmega_WD01.SubGenB.bam
Pmega_WD02.SubGenB.bam
Pmega_WD03.SubGenB.bam
Preni_WH01.SubGenB.bam
Preni_WH02.SubGenB.bam
Preni_WH05.SubGenB.bam
Preni_WA05.SubGenB.bam
Preni_WA06.SubGenB.bam
Preni_WA07.SubGenB.bam
Preni_WA08.SubGenB.bam
Preni_WA09.SubGenB.bam
Preni_WA10.SubGenB.bam
Preni_WA11.SubGenB.bam
Pver_col_WE06.SubGenB.bam
Pver_ver_WE02.SubGenB.bam
Pver_col_WH03.SubGenB.bam
Pver_col_WD09.SubGenB.bam
Pver_macro_WF06.SubGenB.bam
Pver_macro_WE07.SubGenB.bam
Pver_ver_WE12.SubGenB.bam
Pver_macro_WA05.SubGenB.bam
Pver_macro_WH05.SubGenB.bam
Pver_macro_WG06.SubGenB.bam
Pver_macro_WF07.SubGenB.bam
Pver_macro_WH04.SubGenB.bam
Pver_macro_WG05.SubGenB.bam
Pvul_vul_WF05.SubGenB.bam
Pvul_hetero_WA12.SubGenB.bam
Pvul_hetero_WA01.SubGenB.bam
Pvul_hetero_WA02.SubGenB.bam
Pvul_hetero_WA03.SubGenB.bam
Pvul_vul_WG02.SubGenB.bam
Pvul_vul_WG09.SubGenB.bam
Pvul_vul_WG10.SubGenB.bam
Pvul_vul_WF11.SubGenB.bam
Pvul_vul_WE12.SubGenB.bam
Pvul_vul_WE01.SubGenB.bam
Pvul_vul_WE02.SubGenB.bam
Pvul_vul_WE03.SubGenB.bam
Pvul_vul_WE04.SubGenB.bam
Pvul_vul_WE05.SubGenB.bam
Pvul_vul_WE06.SubGenB.bam
Pvul_sib_WC02.SubGenB.bam
Pvul_sib_WC03.SubGenB.bam
Pvul_sib_WC04.SubGenB.bam
Pvul_sib_WC05.SubGenB.bam
Pvul_sib_WC01.SubGenB.bam
Pvul_sib_WB03.SubGenB.bam
grandis_WA05.bam
grandis_WB03.bam
grandis_WB04.bam
grandis_WB10.bam
grandis_WC02.bam
grandis_WC09.bam
grandis_WD08.bam
grandis_WE08.bam
grandis_WF03.bam
grandis_WG03.bam
::::::::::::::
EXCLUDE_FILES/OUTRGROUPS_B.txt
::::::::::::::
Pver_macro_WA05.SubGenB.bam
Pver_macro_WH05.SubGenB.bam
Pver_macro_WG06.SubGenB.bam
Pver_macro_WF07.SubGenB.bam
Pver_macro_WH04.SubGenB.bam
Pver_macro_WG05.SubGenB.bam
::::::::::::::
EXCLUDE_FILES/OUTRGROUPS_B_Pmacro.txt
::::::::::::::
Pvul_vul_WG02.SubGenB.bam
Pvul_vul_WG09.SubGenB.bam
Pvul_vul_WG10.SubGenB.bam
Pvul_vul_WF11.SubGenB.bam
Pvul_vul_WE12.SubGenB.bam
Pvul_vul_WE01.SubGenB.bam
```

```R
#################################################################
### Section Primula (automated)
#################################################################

library(tidyverse)
library(ggpubr)
library(ggplot2)

parse_mk_file <- function(x){
    full_join(x, tsinfo, by="transcript") %>% filter(is_longest=="yes") %>% select(gene, cds_length, pN, pS, dN, dS, f0, f2, f3, f4, pval, dos) %>%
    mutate(logP = -1 * log10(pval), direction = ifelse(dos < 0, "neg", "pos"), sig = ifelse(logP>= 2, "yes", "no"), 
                result = paste0(direction, ":", sig)) %>% # ADDS ONE COLUMN WITH THE logP VALUE AND THREE CONDITIONALS COLUMNS BASED ON THE SIGNIFICANCE 
    filter(!is.na(dos), cds_length <= 10000, !is.na(logP), dN+dS > 2, pN+pS > 2) %>%  #FILTERS SHORT TRANSCRIPTS, NA's and some dN,dS AND pNpS threshold
    mutate(gene_label = str_remove(gene, "gene-")) %>% #CREATES GENE LABEL
    mutate(gene_label = ifelse(str_detect(gene_label, "LOC"), "", gene_label))
}

volcano_plot <- function(x){
    ggplot(data=x, aes(y=logP, x=dos, color=result, alpha=sig)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_manual(values=c("grey", "black",  "grey", "black"), guide="none") +
    scale_alpha_discrete(range = c(0.2, 1), guide="none") +
#    scale_color_manual(values=c("#1f78b4", "#1f78b4",  "#e31a1c", "#e31a1c"), guide="none") +
#    scale_alpha_discrete(range = c(0.2, 0.8), guide="none") +
    xlab("") + 
    scale_x_continuous(limits = c(-1,1)) +
    ylab("-log10 P-value (MK test)") + scale_y_continuous(limits = c(0,3.5))}
#+
#    geom_text(aes(label=ifelse(logP>2,as.character(gene_label),'')), hjust=-0.1,vjust=0.4,color="black") #ADDS LABEL MK-test IS SIGNIFICANT
#}

tsinfo <- read_tsv("transcript-counts.tsv")

mk_Pela_palla <- read_tsv("mk_Pela_palla.tsv")
mk_Pmega <- read_tsv("mk_Pmega.tsv")
mk_Preni <- read_tsv("mk_Preni.tsv")
#mk_Pver_macro <- read_tsv("mk_Pver_macro.tsv")
mk_Pvul_vul <- read_tsv("mk_Pvul_vul.tsv")
mk_Pgrandis_sub_A <-read_tsv("mk_sub_A_veris_5.tsv")
mk_Pgrandis_sub_B <-read_tsv("mk_sub_B_veris_5.tsv")

mk_Pela_palla_filtered <- parse_mk_file(mk_Pela_palla)
mk_Pmega_filtered <- parse_mk_file(mk_Pmega)
mk_Preni_filtered <- parse_mk_file(mk_Preni)
#mk_Pver_macro_filtered <- parse_mk_file(mk_Pver_macro)
mk_Pvul_vul_filtered <- parse_mk_file(mk_Pvul_vul)
mk_Pgrandids_sub_A_filtered <-parse_mk_file(mk_Pgrandis_sub_A)
mk_Pgrandids_sub_B_filtered <-parse_mk_file(mk_Pgrandis_sub_B)

mk_Pela_palla_plot <- volcano_plot(mk_Pela_palla_filtered)
mk_Pmega_plot <- volcano_plot(mk_Pmega_filtered)
mk_Preni_plot <- volcano_plot(mk_Preni_filtered)
#mk_Pver_macro_plot <- volcano_plot(mk_Pver_macro_filtered)
mk_Pvul_vul_plot <- volcano_plot(mk_Pvul_vul_filtered)
mk_Pgrandis_sub_A_plot <- volcano_plot(mk_Pgrandids_sub_A_filtered)
mk_Pgrandis_sub_B_plot <- volcano_plot(mk_Pgrandids_sub_B_filtered)

ggarrange(mk_Pela_palla_plot, mk_Pmega_plot,mk_Preni_plot,
          mk_Pvul_vul_plot,mk_Pgrandis_sub_A_plot,mk_Pgrandis_sub_B_plot,
          labels = c("A","B","C","D","E","F"),
          ncol = 2, nrow = 3)
```
