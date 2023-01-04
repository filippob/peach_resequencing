#!/bin/bash
#########################
## PCA (principal comp.)
#########################

## software
plink=/home/biscarinif/software/plink/plink
admixture=/home/biscarinif/software/admixture_linux-1.3.0/admixture

## input/output paths
homefolder=/home/freeclimb
INPUT_FILE=$homefolder/Analysis/merged_vcf/imputed.vcf.gz
outdir=$homefolder/Analysis/admixture

## parameters
k=20

if [ ! -d ${outdir} ];
then
	## create output dir if it does not exist
	echo "Creating the output folder ${outdir}"
	mkdir -p $outdir
fi


## convert vcf to Plink ped/map
$plink --allow-extra-chr --vcf ${INPUT_FILE} --double-id --recode --out $outdir/temp

for i in {1..8}
do
   echo "Recoding chromosome n. $i"
   sed -i "s/Pp0$i/$i/g" $outdir/temp.map
done

## convert to binary Plink for input to Admixture
$plink --allow-extra-chr --file $outdir/temp --double-id --make-bed --out $outdir/imputed


#############################
## ADMIXTURE		  ###
#############################
echo " - running Admixture ... "

cd $outdir
#$admixture $outdir/imputed.bed $k

## cross validation to find the best value for k (start from 2, 1 doesn't make sense)
for (( K=2; K<=$k; K++ ))  
do 
	echo "trying $K as value for K"
	$admixture --cv $outdir/imputed.bed $K | tee log${K}.out 
done

## clean repo
rm $outdir/temp*

echo "DONE!"
