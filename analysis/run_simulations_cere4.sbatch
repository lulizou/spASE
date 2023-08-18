#!/bin/bash

#SBATCH -N 1 # Number of Nodes
#SBATCH --job-name=run_spase
#SBATCH --output=run_spase_%A_%a.out
#SBATCH --error=run_spase_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --array=1,17,45


module load R
Rscript 02.5_simulation_cere4.R ${SLURM_ARRAY_TASK_ID}

