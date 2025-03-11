#!/bin/bash

#SBATCH --partition=shared
#SBATCH -J abrEgg
#SBATCH --output Abrown_egg.log
#SBATCH --open-mode=append
##SBATCH --mem=512000
module load anaconda/colsa
conda activate eggnog-mapper-2.1.12 

for assembly in megahit_*/final.contigs.fa
do
sp=${assembly%/*}
emapper.py -m diamond --cpu 24 --itype metagenome --translate --genepred prodigal  -i $assembly -o $sp
