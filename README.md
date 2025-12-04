# Myosins generate contractile force and maintain organization in the cytokinetic contractile ring
Authors: Zachary McDargh, Shuyuan Wang, Harvey F. Chin, Sathish Thiyagarajan,  View ORCID ProfileErdem Karatekin, Thomas D. Pollard, Ben O’Shaughnessy
doi: https://doi.org/10.1101/2021.05.02.442363

## Running the code
### Option 1: running the code in an HPC cluster
To run the code in an HPC cluster, please do "sbatch runv223.sh".

### Option 2: running the code in a local computer
To run the code in a local computer, please open 'task_yeti_v223.m' using MATLAB_R2024b. At the beginning of the Code, please define itask to be any random number. 

## Result file
Results will be saved as a series of compiled .mat files depending on t. t is the constriction time. The files of results will be in the same location as 'task_yeti_v233.m'. The name of the files will be 'wtcon_{itask}_{t}sec.mat'.

## Analysis
The code used to analyze the data is located in the analysis repository. The filenames correspond to their respective functions.
