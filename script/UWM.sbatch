#!/bin/sh -e

#SBATCH --array=1-7
#SBATCH --output=case-%a.out
#SBATCH --error=case-%a.err

points=1000000000

# Make sure last argument matches total number of jobs in array!

./PoreDy < input-$SLURM_ARRAY_TASK_ID > Case$SLURM_ARRAY_TASK_ID.TD
