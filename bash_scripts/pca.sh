#!/bin/bash
#########################
## PCA (principal comp.)
#########################

## command lines (and instruction) for data filtering, imputation and preparation (for GWAS)
plink=/home/filippo/Downloads/plink

homefolder=$HOME/Documents/freeclimb/VariantCalling
INPUT_FILE=$homefolder/Analysis/imputed.vcf.gz
fname="$(basename -- $INPUT_FILE)"
outdir=$homefolder/Analysis

#########################
## PCA
#########################

## imputation of missing genotypes

## convert to vcf
$plink --allow-extra-chr --vcf ${INPUT_FILE} --double-id --pca --out $outdir/$(basename -- $INPUT_FILE)

## clean repo
rm $outdir/*.log $outdir/*.nosex
