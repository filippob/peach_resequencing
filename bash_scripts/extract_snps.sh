#!/bin/bash
set -euxo pipefail

## script to read an input vcf file with SNP only (output of the snp_only bash script) and a list of SNP (chromosome, position),
## and then extract the SNPs in the list (if present) from the input vcf file
## this script can be used before or after filtering/merging input vcf file(s)

### files and folders
input_file=$1 ## e.g. Analysis/IGA/results/snps_only.vcf.gz
seqrun=$3 ## e.g. IGA
snp_list=$2 ## e.g. /home/freeclimb/Analysis/merged_vcf/lista_SNP.txt
homefolder="/home/freeclimb"
#dir="$(dirname "${input_file}")"
#outdir="Analysis/${seqrun}/filtered"
outdir="Analysis/$seqrun/extracted_vcf"

### software
bcftools_img='/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-bcftools:1.15.1.img'
vcftools_img='/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-vcftools0.1.16.img'

### parameters

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
temp="$(basename "${input_file}")"
fname=${temp/.vcf.gz/""}
outfile="$homefolder/$outdir/${fname}_filtered"
echo "will be writing out to file $outfile"

if [ -f "$outfile" ]; then
    	
	echo "$outfile exists. Skip the filtering step"

else
	echo " - extracting SNPs in the provided list form the input Vcf file"
	#singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools view -e 'QUAL <= 1 || INFO/DP > 35 || INFO/DP < 0.2' -c $mac -q $maf ${input_file} | bgzip > $outfile
	singularity run -B /home/freeclimb/:/home/freeclimb/ ${vcftools_img} vcftools --gzvcf ${input_file} --positions ${snp_list} --recode --recode-INFO-all --out $outfile
	echo "compressing vcf file wth bgzip"
	singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bgzip -c $outfile.recode.vcf > $outfile.vcf.gz
	echo " - indexing the subset vcf file with tabix"
	singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} tabix ${outfile}.vcf.gz
	rm ${outfile}.recode.vcf
fi

## counting n. of samples
nsamples=`singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools query -l $outfile.vcf.gz | wc -l`
echo "n. of output samples from $outfile is $nsamples"

## counting n. of variants
nvars=`singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools query -f '%POS\n' $outfile.vcf.gz | wc -l`
echo "n. of output variants from $outfile is $nvars"

## stat report from bcftools
echo " - calculating stats on filtered file ... "
singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools stats $outfile.vcf.gz > ${outdir}/bcftools.filtered.stats

echo "DONE!!"
