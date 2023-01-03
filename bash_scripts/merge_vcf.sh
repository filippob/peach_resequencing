#!/bin/bash
set -euxo pipefail

## script to read an input vcf file with SNP only (output of the snp_only bash script) and
## filter SNP sites according to multiple criteria (see below); results are written to a new vcf file and tabix-indexed

### files and folders
#input_file=$1 ## e.g. Analysis/IGA/results/snps_only.vcf.gz
homefolder="/home/freeclimb"
f1="$homefolder/Analysis/IGA/filtered/snps_only_filtered.vcf.gz"
f2="$homefolder/Analysis/BGI/filtered/snps_only_filtered.vcf.gz"
f3="$homefolder/Analysis/Reseq/filtered/snps_only_filtered.vcf.gz"
f4="$homefolder/Analysis/NCBI/filtered/snps_only_filtered.vcf.gz"
#dir="$(dirname "${input_file}")"
outdir="Analysis/merged_vcf"

### software
bcftools_img='/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-bcftools:1.15.1.img'
vcftools_img='/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-vcftools0.1.16.img'

### parameters


echo "-----------------------"
echo "MERGING VCF FILES"
echo "-----------------------"

echo "working on the $outdir data experiment"

if [ ! -d $homefolder/${outdir} ];
then
	## create output dir if it does not exist
	echo "Creating the output folder ${outdir}"
	mkdir -p $homefolder/$outdir
fi


echo "bcftools merge"
outfile="$homefolder/$outdir/merged.vcf"
echo "will be writing out to file $outfile"

if [ -f "$outfile" ]; then
    	
	echo "$outfile exists. Skip the filtering step"

else
	echo " - merging multiple  input vcf files"
	singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools merge -m none $f1 $f2 -O v -o ${outfile}
	#singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools view -e 'QUAL <= 1 || INFO/DP > 35 || INFO/DP < 0.2' -c $mac -q $maf ${input_file} | bgzip > $outfile
	echo "compressing vcf file wth bgzip"
	singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bgzip -c $outfile > $outfile.gz
	echo " - indexing the subset vcf file with tabix"
	singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} tabix ${outfile}.gz
fi

## counting n. of samples
nsamples=`singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools query -l $outfile.gz | wc -l`
echo "n. of output samples from $outfile is $nsamples"

## counting n. of variants
nvars=`singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools query -f '%POS\n' $outfile.gz | wc -l`
echo "n. of output variants from $outfile is $nvars"

## stat report from bcftools
echo " - calculating stats on filtered file ... "
singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools stats $outfile.gz > ${outdir}/bcftools.merged.stats

echo "DONE!!"
