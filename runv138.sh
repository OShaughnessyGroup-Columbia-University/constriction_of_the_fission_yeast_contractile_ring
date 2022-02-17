#!/bin/sh

# Directives
#PBS -N transient
#PBS -W group_list=yeticheme
#PBS -l nodes=1:ppn=1,walltime=72:00:00,mem=8gb
#PBS -M zm2296@columbia.edu
#PBS -m n
#PBS -V
#PBS -t 1-24

# Set output and error directories
#PBS -o localhost:/vega/cheme/users/zm2296/outputs
#PBS -e localhost:/vega/cheme/users/zm2296/outputs

matlab -nosplash -nodisplay -nodesktop -r "task_yeti_v138($PBS_ARRAYID)" 

# End of script
