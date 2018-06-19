# RNA-seq analysis of R178 data set

For RNA-seq analysis [RNAcocktail](https://bioinform.github.io/rnacocktail/) pipelines were used.

Pipelines:

1. raw_reads -> SALMON-SMEM -> DESeq2
2. raw_reads -> HISAT2 -> Stringtie -> Ballgown

For first pipeline two results are obtained: with and with out quality check control for reads. Pipeline #1 with out quality check was performed on ensembl transcriptome annotation version 67 (version 92 is available). 