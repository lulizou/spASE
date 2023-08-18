#!/bin/bash

#SBATCH -N 1 # Number of Nodes
#SBATCH -c 8 # Number of threads
#SBATCH --job-name=run_spasev1_combined
#SBATCH --output=run_spasev1_combined_%A_%a.out
#SBATCH --error=run_spasev1_combined_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --array=4

ARRAY=('run_spasev1_combined_hippo1.R' 'run_spasev1_combined_hippo3.R' 'run_spasev1_combined_cere4.R' 'run_spase_hippo2.R' 'run_spasev1_visium_cere.R' 'run_spase_mixture.R')

module load R
Rscript ${ARRAY[${SLURM_ARRAY_TASK_ID}]}

