#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH -A cheme                # The account name for the job.
#SBATCH -J clean                # The job name.
#SBATCH -t 1:00              # The time the job will take to run.
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb       # The memory the job will use per cpu core.
#SBATCH --array=3,7,8,9,10,15,19,21,22,23,27,28,30,32,33,35,41,43,46,48,51,54,55,56,60
rm -fv brx6*_$SLURM_ARRAY_TASK_ID*.mat

# End of script
