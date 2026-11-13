#!/bin/bash
#SBATCH --job-name=conduit_timing
#SBATCH --account=jkarvane
#SBATCH --partition=longrun
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-100
#SBATCH --output=logs/timing_%A_%a.out
#SBATCH --error=logs/timing_%A_%a.err
 
 
module load r-env/461
 
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
 
srun Rscript --no-save 07_polynomial_simulation.R $SLURM_ARRAY_TASK_ID
