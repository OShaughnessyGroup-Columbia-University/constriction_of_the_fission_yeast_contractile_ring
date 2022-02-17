#!/bin/sh

# Directives
#PBS -N wtfc3_zm2
#PBS -W group_list=yeticheme
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=4gb
#PBS -M zm2296@columbia.edu
#PBS -m n
#PBS -V
#PBS -t 1-5

# Set output and error directories
#PBS -o localhost:/vega/cheme/users/zm2296/outputs
#PBS -e localhost:/vega/cheme/users/zm2296/outputs

# detached ring simulation

matlab -nosplash -nodisplay -nodesktop -r "task_yeti_v111_zm2($PBS_ARRAYID)" 

# End of script
