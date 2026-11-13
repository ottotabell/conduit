#!/bin/bash
#SBATCH --job-name=causal_sim
#SBATCH --account=jkarvane
#SBATCH --partition=longrun
#SBATCH --time=200:00:00
#SBATCH --array=1-18                # = length(networks)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --output=logs00/sim_%A_%a.out
#SBATCH --error=logs00/sim_%A_%a.err

module load r-env

srun Rscript --no-save 05_simulation00.R $SLURM_ARRAY_TASK_ID