#!/bin/bash

#SBATCH -N 1 # Number of Nodes
#SBATCH -c 8 # Number of threads
#SBATCH --job-name=run_spasev1_combined
#SBATCH --output=run_spasev1_combined_%A_%a.out
#SBATCH --error=run_spasev1_combined_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --array=0,1,4

ARRAY=('scripts/run_spase_hippo1_overall_spatial.R' 'scripts/run_spase_hippo2_overall_spatial.R' 'scripts/run_spase_hippo3_overall_spatial.R' 'scripts/run_spase_cere3_overall_spatial.R' 'scripts/run_spase_cere4_visium_overall_spatial.R' 'scripts/run_spase_mix5_visium_overall_spatial.R')
module load R
Rscript ${ARRAY[${SLURM_ARRAY_TASK_ID}]}

