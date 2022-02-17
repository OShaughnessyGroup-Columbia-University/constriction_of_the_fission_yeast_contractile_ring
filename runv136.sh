#!/bin/sh

# Directives
#PBS -N bind2d8
#PBS -W group_list=yeticheme
#PBS -l nodes=1:ppn=1,walltime=4:00:00,mem=4gb
#PBS -M zm2296@columbia.edu
#PBS -m n
#PBS -V
#PBS -t 8

# Set output and error directories
#PBS -o localhost:/vega/cheme/users/zm2296/outputs
#PBS -e localhost:/vega/cheme/users/zm2296/outputs

# detached ring simulation

matlab -nosplash -nodisplay -nodesktop -r "task_yeti_v136($PBS_ARRAYID)" 

# End of script
