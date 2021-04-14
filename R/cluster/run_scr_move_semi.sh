#!/bin/bash
#SBATCH --partition=cymerman_p
#SBATCH --job-name=scr_move_semi
#SBATCH --time=10-01:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=rbchan@uga.edu
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

module load R/4.0.0-foss-2019b

time R CMD BATCH --vanilla scr_move_semi.R scr_move_semi.Rout
