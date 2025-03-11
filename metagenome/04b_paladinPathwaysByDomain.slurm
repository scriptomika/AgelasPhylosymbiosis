#!/bin/bash

#SBATCH --partition=shared
#SBATCH -J pcountA
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1 
#SBATCH --output pcountsA.log
#SBATCH -x node117,node118

module load linuxbrew/colsa

# SYNOPSIS
# 1. run initial @@pathways command to identify all pathways present
# 2. pull raw counts mapping to those pathways, with respect to Bacterial, Archaeal or Metazoan origin

#outputs 3 files per pathway per library

#pull all ec pathway accessions identified from across samples:
#only run once all initial @@pathways runs completed
# cat post_pathways/*_pathways.txt |cut -f1 |sort |uniq | head -n-1 > ec_accs.txt

for tsv in post_aggreg/*tsv 
do 

libpath=${tsv%%.both_uniprot.tsv}
lib=${libpath#post_aggreg/}
echo $lib

echo "Starting counts..."; startGO=$(date  +%T); echo "$startGO"

#loop thru all the ec accessions for Bacteria, Archaea, Metazoa

while IFS= read -r ecacc; do
if [ ! -f post_pathways/counts/${lib}_Bacteria_${ecacc}_counts.txt ]
then
echo "Getting counts across Bacteria for pathway $ecacc in sample $lib"

HTTP_PROXY=HTTPS_PROXY=http://premise.sr.unh.edu:3128 paladin-plugins.py @@pathways \
-i $tsv -q 20 -p $ecacc -s Bacteria  @@write ${lib}_Bacteria_${ecacc}_counts.txt
mv ${lib}*_counts.txt post_pathways/counts/
fi

if [ ! -f post_pathways/counts/${lib}_Archaea_${ecacc}_counts.txt ]
then
echo "Getting counts across Archaea for pathway $ecacc in sample $lib"

HTTP_PROXY=HTTPS_PROXY=http://premise.sr.unh.edu:3128 paladin-plugins.py @@pathways \
-i $tsv -q 20 -p $ecacc -s Archaea  @@write ${lib}_Archaea_${ecacc}_counts.txt
mv ${lib}*_counts.txt post_pathways/counts/

fi

if [ ! -f post_pathways/counts/${lib}_Metazoa_${ecacc}_counts.txt ]
then
echo "Getting counts across Metazoa for pathway $ecacc in sample $lib"

HTTP_PROXY=HTTPS_PROXY=http://premise.sr.unh.edu:3128 paladin-plugins.py @@pathways \
-i $tsv -q 20 -p $ecacc -s Metazoa  @@write ${lib}_Metazoa_${ecacc}_counts.txt

mv ${lib}*_counts.txt post_pathways/counts/
fi

done <  <(cat ec_accs.txt )

echo "Finished counting enxymes for $lib."; endGO=$(date  +%T); echo "$endGO"

done
