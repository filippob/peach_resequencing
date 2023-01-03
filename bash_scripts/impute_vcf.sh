#!/bin/bash
set -euxo pipefail

## script to read an input vcf file with SNP only (output of the snp_only bash script) and
## filter SNP sites according to multiple criteria (see below); results are written to a new vcf file and tabix-indexed

### files and folders
input_file=$1 ## e.g. Analysis/IGA/results/snps_only.vcf.gz
homefolder="/home/freeclimb"
#dir="$(dirname "${input_file}")"
outdir="Analysis/merged_vcf"

### software
bcftools_img='/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-bcftools:1.15.1.img'
vcftools_img='/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-vcftools0.1.16.img'
beagle="/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-beagle5.2.img"
plink="/home/biscarinif/software/plink/plink"

### parameters


echo "-------------------------"
echo "IMPUTING MISSING SNP DATA"
echo "-------------------------"

echo "working on the $outdir data experiment"

if [ ! -d $homefolder/${outdir} ];
then
	## create output dir if it does not exist
	echo "Creating the output folder ${outdir}"
	mkdir -p $homefolder/$outdir
fi


echo "Imputation of missing genotypes"
outfile="$homefolder/$outdir/imputed"
echo "will be writing out to file $outfile"

if [ -f "$outfile" ]; then
    	
	echo "$outfile exists. Skip the filtering step"

else
	echo " - using Plink to convert back and to plink/vcf to get rid of contradictory header/format info"
	$plink --allow-extra-chr --vcf $homefolder/${input_file} --chr Pp01,Pp02,Pp03,Pp04,Pp05,Pp06,Pp07,Pp08 --make-bed --out $outdir/temp
        $plink --allow-extra-chr --bfile $outdir/temp --recode vcf-iid --out $outdir/input
	echo " - imputing missing SNP genotypes from vcf file"
	singularity run -B /home/freeclimb/:/home/freeclimb/ $beagle beagle gt=$outdir/input.vcf out=$outfile
	echo "compressing vcf file wth bgzip"
	singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bgzip -c $outfile.vcf.gz > $outfile.vcf.gz
	echo " - indexing the subset vcf file with tabix"
	singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} tabix ${outfile}.vcf.gz
	## janitor
	rm $outdir/temp*
	rm $outdir/input*
fi

## counting n. of samples
nsamples=`singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools query -l $outfile.vcf.gz | wc -l`
echo "n. of output samples from $outfile is $nsamples"

## counting n. of variants
nvars=`singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools query -f '%POS\n' $outfile.vcf.gz | wc -l`
echo "n. of output variants from $outfile is $nvars"

## stat report from bcftools
echo " - calculating stats on filtered file ... "
singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools stats $outfile.vcf.gz > ${outdir}/bcftools.imputed.stats

echo "DONE!!"
