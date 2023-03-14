#!/bin/bash
#SBATCH --job-name=GxE_markers   		# Job name
#SBATCH --get-user-env
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single task
#SBATCH --cpus-per-task=1             # Declare 1 CPUs per task
#SBATCH --mem=4gb                     # Job memory request
#SBATCH --output=/home/freeclimb/log/GxE_markers_MD_2020_%j.log   # Standard output and error log

## !! REMEMBER: CHANGE CONFIG FILE TO USE THE BEAGLE SINGULARITY CONTAINER !!
repofolder="/home/freeclimb/peach_resequencing"
start=`date +%s`
Rscript $repofolder/r_scripts/GXE_bglr_case1_markers.R
end=`date +%s`
runtime=$((end-start))

echo "execution time was $runtime"
echo "DONE!"
