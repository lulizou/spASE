#!/bin/bash

#SBATCH -N 1 # Number of Nodes
#SBATCH --job-name=run_spase
#SBATCH --output=run_spase_%A_%a.out
#SBATCH --error=run_spase_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH --array=3-18


module load R
Rscript 02.5_simulation_cere3.R ${SLURM_ARRAY_TASK_ID}

