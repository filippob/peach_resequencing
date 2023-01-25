# peach_resequencing
Scripts for variant calling from resequencing data on peach accessions (Prunus persica)

- make sample sheet for IGA samples: R script `make_samplesheet.r`
- make sample sheet for BGI samples: Python script `make_sample_sheet.py` (to be completed: include the case of multiple files within subfolders --> more than two files per samples, multiple R1/R2 fastq files)

Bash scripts to use `bcftools` and `vcftools` to subset the vcf file from variant calling and to produce summary statistics

1. normalise_vcf.sh: sometimes Freebayes adds the whole haplotype around a SNP --> need to remove the non-polymorphic sites and leave only yhe SNP
2. snp_only_from_vcf.sh: take only biallelic SNPs (remove MNPs/Indels, triallelic sites)
3. stats: frequency, missing rate
4. quality filter: how?
5. merge?

## GxE interactions

1. approaches with BGLR from https://github.com/gdlc/BGLR-R/blob/master/inst/md/GxE_usingInteractions.md (paper: https://pubmed.ncbi.nlm.nih.gov/25660166/)
