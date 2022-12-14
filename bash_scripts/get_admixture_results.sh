#!/bin/sh

## script to collect results from cross-validation with Admixture

## input/output paths
homefolder=/home/filippo/Documents/freeclimb/VariantCalling
outdir=$homefolder/Analysis/admixture
k=20

cd $outdir
#get the cross-validation error (-h avoids the file name to be printed each time)
grep -h "CV" log*.out > CVerr

#get the number of iterations to reach convergence (start from 2, 1 cluster makes no sense - unless it was run by mistake in admixture.sh)
for (( i=2; i<=$k; i++ )); do cat log$i.out | awk '{if (NF==8 && $2=="(QN/Block)") print $1}' | tail -1; done > niter

echo "DONE!"

