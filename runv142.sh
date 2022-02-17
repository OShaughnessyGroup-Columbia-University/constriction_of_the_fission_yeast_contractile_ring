#!/bin/sh

# Directives
#PBS -N weak_myp_unbind
#PBS -W group_list=yeticheme
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=4gb
#PBS -M zm2296@columbia.edu
#PBS -m n
#PBS -V
#PBS -t 1-10

# Set output and error directories
#PBS -o localhost:/vega/cheme/users/zm2296/outputs
#PBS -e localhost:/vega/cheme/users/zm2296/outputs

matlab -nosplash -nodisplay -nodesktop -r "task_yeti_v142($PBS_ARRAYID)" 

# End of script