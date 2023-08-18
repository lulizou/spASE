#!/bin/bash

#SBATCH -N 1 # Number of Nodes
#SBATCH -c 8 # Number of threads
#SBATCH --job-name=run_spase_combined
#SBATCH --output=run_spase_combined_%A_%a.out
#SBATCH --error=run_spase_combined_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --array=4

ARRAY=('run_spase_combined_hippo1.R' 'run_spase_combined_hippo3.R' 'run_spase_combined_cere4.R' 'run_spase_hippo2.R' 'run_spase_cere4_visium.R')

module load R
Rscript ${ARRAY[${SLURM_ARRAY_TASK_ID}]}

