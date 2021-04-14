#!/bin/bash
#SBATCH --partition=cymerman_p
#SBATCH --job-name=session_info
#SBATCH --time=01:00:00
#SBATCH --mem=1G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

cd $SLURM_SUBMIT_DIR

module load R/4.0.0-foss-2019b

time R CMD BATCH --vanilla session_info.R
