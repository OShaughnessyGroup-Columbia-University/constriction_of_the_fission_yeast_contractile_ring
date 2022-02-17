#!/bin/sh

# Directives
#PBS -N task1
#PBS -W group_list=yeticheme
#PBS -l nodes=1:ppn=1,walltime=72:00:00,mem=8gb
#PBS -M sw2703@columbia.edu
#PBS -m n
#PBS -V
#PBS -t 1-10

# Set output and error directories
#PBS -o localhost:/vega/cheme/users/sw2703/outputs
#PBS -e localhost:/vega/cheme/users/sw2703/outputs

# detached ring simulation

matlab -nosplash -nodisplay -nodesktop -r "task_yeti_v32($PBS_ARRAYID)" 

# End of script