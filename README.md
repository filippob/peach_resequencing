# peach_resequencing
Scripts for variant calling from resequencing data on peach accessions (Prunus persica)

- make sample sheet for IGA samples: R script `make_samplesheet.r`
- make sample sheet for BGI samples: Python script `make_sample_sheet.py` (to be completed: include the case of multiple files within subfolders --> more than two files per samples, multiple R1/R2 fastq files)

Bash scripts to use `bcftools` and `vcftools` to subset the vcf file from variant calling and to produce summary statistics

1. normalise_vcf.sh: sometimes Freebayes adds the whole haplotype around a SNP --> need to remove the non-polymorphic sites and leave only yhe SNP
2. filter (take only SNP)
