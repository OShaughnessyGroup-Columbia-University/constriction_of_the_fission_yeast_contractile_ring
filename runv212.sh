#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH -A cheme                # The account name for the job.
#SBATCH -J wrax                  # The job name.
#SBATCH -t 6:00:00              # The time the job will take to run.
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2gb       # The memory the job will use per cpu core.
#SBATCH --array=1-132

module load matlab/2018b
matlab -nosplash -nodisplay -nodesktop -r "task_yeti_v212($SLURM_ARRAY_TASK_ID)" 

# End of script

