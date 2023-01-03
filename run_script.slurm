#!/bin/bash
#SBATCH --job-name=filter-vcf   		# Job name
#SBATCH --get-user-env
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single task
#SBATCH --cpus-per-task=4             # Declare 1 CPUs per task
#SBATCH --mem=24gb                     # Job memory request
#SBATCH --output=/home/freeclimb/log/filter_vcf_%j.log   # Standard output and error log

## PARAMETERS
home="/home/freeclimb"
repofolder="${home}/peach_resequencing"
inpfile="${home}/Analysis/BGI/snp_only/snps_only.vcf.gz"
label="BGI" ## used for subfolder

## RUN
start=`date +%s`
bash $repofolder/bash_scripts/filter_vcf.sh $inpfile $label
end=`date +%s`
runtime=$((end-start))

echo "execution time was $runtime"
echo "DONE!"

