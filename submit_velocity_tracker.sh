#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH -A cheme                # The account name for the job.
#SBATCH -J trck                 # The job name.
#SBATCH -t 2:00:00              # The time the job will take to run.
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4gb       # The memory the job will use per cpu core.
#SBATCH --array=0-11

module load matlab/2018b
INDX=(1 11 21 31 41 51 2 12 22 32 42 52)
TINDX=(1200:3600 1200:3600 1200:3600 300:900 300:900 300:900 1200:3600 1200:3600 1200:3600 300:900 300:900 300:900)
matlab -nosplash -nodisplay -nodesktop -r "node_velocity_tracker(${INDX[$SLURM_ARRAY_TASK_ID]}, 'coopf','dump',${TINDX[$SLURM_ARRAY_TASK_ID]})" 

# End of script
