#!/bin/bash

#SBATCH -N 1 # Number of Nodes
#SBATCH -c 8 # Number of threads
#SBATCH --job-name=run_betabin
#SBATCH --output=run_betabin_%A_%a.out
#SBATCH --error=run_betabin_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --array=0-1

ARRAY=('run_betabin_bias_and_celltype_cere_4.R' 'run_betabin_bias_and_celltype_hippo_3.R')
module load R
Rscript ${ARRAY[${SLURM_ARRAY_TASK_ID}]}

