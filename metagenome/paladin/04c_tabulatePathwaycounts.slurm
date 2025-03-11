#!/bin/bash

#SBATCH --partition=shared
#SBATCH -J pathwy
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1 
#SBATCH --output msp_tabulate.log
#SBATCH -x node117,node118

module load linuxbrew/colsa
Rscript --vanilla ./04c_tabulate_counts.R
