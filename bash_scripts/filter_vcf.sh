#!/bin/bash
set -euxo pipefail

## script to read an input vcf file from the Nextflow resequencing-mem pipeline and
## keep only BIALLELIC SNPs; results are written to a new vcf file and tabix-indexed

### files and folders
input_file=$1
seqrun=$2
homefolder="/home/freeclimb"
#dir="$(dirname "${input_file}")"
outdir="Analysis/${seqrun}/filtered"

### software
bcftools_img='/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-bcftools:1.15.1.img'
vcftools_img='/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-vcftools0.1.16.img'

### parameters
thin=1e+6 ## distance between adjacent sites
maf=0.05 ## 0.05 --> 5%
mac=8 ## min. number of allele counts
minDP=0.2 ## average coverage across all samples
maxmiss=0.5 ## max missing rate per site


echo "-----------------------"
echo "FILTERING SNP VCF FILES"
echo "-----------------------"

echo "working on the $outdir data experiment"

if [ ! -d $homefolder/${outdir} ];
then
	## create output dir if it does not exist
	echo "Creating the output folder ${outdir}"
	mkdir -p $homefolder/$outdir
fi


## stat report from bcftools on the input file
singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools stats ${input_file} > ${outdir}/bcftools.input.stats

echo "analysing data from ${input_file}"
fname="$(basename "${input_file}")"
outfile="$homefolder/$outdir/${fname}_filtered.vcf.gz"
echo "will be writing out to file $outfile"

if [ -f "$outfile" ]; then
    	
	echo "$outfile exists. Skip the filtering step"

else
	echo " - extracting only biallelic SNPs from the input vcf file"
	singularity run -B /home/freeclimb/:/home/freeclimb/ ${vcftools_img} vcftools --gzvcf ${input_file} --remove-indels --thin $thin --maf $maf --mac $mac --max-missing $maxmiss --min-meanDP $minDP --recode --recode-INFO-all --out $outfile

	echo " - indexing the subset vcf file with tabix"
	singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} tabix $outfile
fi

## counting n. of samples
nsamples=`singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools query -l $outfile | wc -l`
echo "n. of output samples from $outfile is $nsamples"

## counting n. of variants
nvars=`singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools query -f '%POS\n' $outfile | wc -l`
echo "n. of output variants from $outfile is $nvars"

## stat report from bcftools
singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools stats $outfile > ${outdir}/bcftools.filtered.stats

echo "DONE!!"
