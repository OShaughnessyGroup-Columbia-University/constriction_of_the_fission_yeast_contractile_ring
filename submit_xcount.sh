#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH -A cheme                # The account name for the job.
#SBATCH -J xcnt                 # The job name.
#SBATCH -t 10:00                # The time the job will take to run.
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb       # The memory the job will use per cpu core.
#SBATCH --array=1-300

module load matlab/2018b
matlab -nosplash -nodisplay -nodesktop -r "count_crossings('rkdt1', $SLURM_ARRAY_TASK_ID, 120)" 

# End of script
