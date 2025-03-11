#!/bin/bash

#SBATCH --partition=shared
#SBATCH -J agg
#SBATCH --output agg.log
#SBATCH --mem 500Gb
#SBATCH --open-mode=append
#SBATCH -x node117,node118
module load linuxbrew/colsa


# aggregates _1 and _2 paladin results for all files in dir 'paladin_results'

for sam in paladin_results/*1_paladin.sam
do 
libpath=${sam%_1_paladin.sam}
lib=${libpath#paladin_results/}
samfile=${sam#paladin_results/}

if [ ! -f post_aggreg/${lib}.both_uniprot.tsv ]
then
#echo "no file for $lib"

# aggregate _1 and _2 outputs
paladin-plugins.py @@aggregation -r paladin_results/ -s $samfile -o ${lib}.both
mv ${lib}.both* post_aggreg/
fi
done
