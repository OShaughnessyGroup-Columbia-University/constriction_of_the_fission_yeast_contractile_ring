#!/bin/sh

# Directives
#PBS -N bind2d5
#PBS -W group_list=yeticheme
#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=8gb
#PBS -M sw2703@columbia.edu
#PBS -m n
#PBS -V
#PBS -t 1-64

# Set output and error directories
#PBS -o localhost:/vega/cheme/users/sw2703/outputs
#PBS -e localhost:/vega/cheme/users/sw2703/outputs

# detached ring simulation

matlab -nosplash -nodisplay -nodesktop -r "task_yeti_v58($PBS_ARRAYID)" 

# End of script