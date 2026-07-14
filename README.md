# *Primula grandis*
This github repository contains the bioinformatic scripts to perform all the population genetic analyses and plots for [Mora-Carrera *et al* (2025)](https://doi.org/10.1093/molbev/msaf162). Genomic patterns of loss of distyly and polyploidization in primroses.

For script regarding the genome assembly, genome scaffolding, gene annotation and transposable elements annotation, please contact [Dr. Narjes Yousefi](mailto:narjes.yousefi2@uzh.ch).

# Preprocessing of population genetic analyses

1) Trimming of short-read sequences with [Trimmomatic](https://github.com/usadellab/Trimmomatic). [Here](/trimmomatic.md).
2) Mapping to reference genome with [BWA](https://github.com/lh3/BWA). [Here](/bwa.md)
3) Variant calling with [BCFtools](https://github.com/samtools/bcftools). [Here](/bcftools.md).

## Genomic analyses results

1) [Subgenome characterization](subphaser.md) (Figure 2b).
2) 
All necessary file to perfom downstream analyses including the reference genome, genome annotation and filtered VCF files can be found in [figshare](https://figshare.com/articles/journal_contribution/Evolutionary_history_of_i_Primula_grandis_i_/28540910)
