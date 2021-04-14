#!/bin/bash
#SBATCH --partition=cymerman_p
#SBATCH --job-name=scr_DA
#SBATCH --time=12-01:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=rbchan@uga.edu
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

module load R/4.0.0-foss-2019b

time R CMD BATCH --vanilla scr_DA.R
