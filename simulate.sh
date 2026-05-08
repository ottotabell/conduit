#!/bin/bash
#SBATCH --job-name=causal_sim
#SBATCH --account=jkarvane
#SBATCH --partition=small          # 'test' for short jobs, 'small' for longer
#SBATCH --time=04:00:00
#SBATCH --array=1-17                # One job per graph — set to length(networks)
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

module load r-env

if test -f ~/.Renviron; then
  sed -i '/TMPDIR/d' ~/.Renviron
fi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs results

srun Rscript simulation.R $SLURM_ARRAY_TASK_ID