#!/bin/bash

#SBATCH -N 1 # Number of Nodes
#SBATCH -c 8 # Number of threads
#SBATCH --job-name=run_rctd
#SBATCH --output=run_rctd_%A_%a.out
#SBATCH --error=run_rctd_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --array=1-5

ARRAY=('./scripts/run_rctd_hippo1.R' './scripts/run_rctd_hippo3.R' './scripts/run_rctd_cere3.R' './scripts/run_rctd_hippo2.R' './scripts/run_rctd_cere4_visium.R' './scripts/run_rctd_mix5_visium.R')

module load R
Rscript ${ARRAY[${SLURM_ARRAY_TASK_ID}]}

