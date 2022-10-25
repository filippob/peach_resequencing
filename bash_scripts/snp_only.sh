#|/bin/bash
set -euxo pipefail

## script to read an input vcf file from the Nextflow resequencing-mem pipeline and
## keep only BIALLELIC SNPs

### files and folders
input_file=$1
homefolder="/home/freeclimb"
outdir="Analysis/IGA/results"

### software
bcftools_img='/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-bcftools:1.15.1.img'

### parameters
vartype="snps"
min_allele=2 ## min n. of alleles at SNP site
max_allele=2 ## max n. of alleles at SNP site


echo "----------------------"
echo "KEEPING BIALLELIC SNPs"
echo "----------------------"

echo "working on the $outdir data experiment"

if [ ! -d $homefolder/${outdir} ];
then
	## create output dir if it does not exist
	echo "Creating the output folder ${outdir}"
	mkdir -p $homefolder/$outdir
fi


echo "analysing data from ${input_file}"
outfile="$homefolder/$outdir/snps_only.vcf.gz"
echo "will be writing out to file $outfile"

echo " - extracting only biallelic SNPs from the input vcf file"
singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} bcftools view --types $vartype -m ${min_allele} -M ${max_allele} ${input_file} --output $outfile

echo " - indexing the subset vcf file with tabix"
singularity run -B /home/freeclimb/:/home/freeclimb/ ${bcftools_img} tabix $outfile

echo "DONE!!"
