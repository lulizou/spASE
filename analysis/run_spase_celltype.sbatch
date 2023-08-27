#!/bin/bash

#SBATCH -N 1 # Number of Nodes
#SBATCH -c 8 # Number of threads
#SBATCH --job-name=run_spase_combined
#SBATCH --output=run_spase_combined_%A_%a.out
#SBATCH --error=run_spase_combined_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --array=0-4

ARRAY=('scripts/run_spase_hippo1_celltype.R' 'scripts/run_spase_hippo2_celltype.R' 'scripts/run_spase_hippo3_celltype.R' 'scripts/run_spase_cere3_celltype.R' 'scripts/run_spase_cere4_visium_celltype.R')
module load R
Rscript ${ARRAY[${SLURM_ARRAY_TASK_ID}]}

