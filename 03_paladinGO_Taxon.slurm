#!/bin/bash

#SBATCH --partition=shared
#SBATCH -J agpostGO
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1 
#SBATCH --output agelpostGO.log
#SBATCH -x node117,node118

module load linuxbrew/colsa


# AGELAS

#for tsv in paladin_results/*_1_paladin_uniprot.tsv
for tsv in post_aggreg/*.both_uniprot.tsv

do 

libpath=${tsv%%.both_uniprot.tsv}
lib=${libpath#post_aggreg/}
echo $lib
if [ ! -f post_GO/${lib}* ]
then

echo "Starting GO step..."; startGO=$(date  +%T); echo "$startGO"
#group GO terms and write abundances, filtering for mapping quality >= 20
HTTP_PROXY=HTTPS_PROXY=premise.sr.unh.edu:3128 paladin-plugins.py @@go -i $tsv -q 20 @@write ${lib}_GO_output.txt
#echo "Finished GO step."; endGO=$(date  +%T); echo "$endGO"
mv ${lib}_GO_output.txt post_GO/

# get Bacteria taxomomic breakdown
echo "Starting Bacterial taxonomy step..."; starttax1=$(date +%T); echo "$starttax1"
paladin-plugins.py @@taxonomy -i $tsv -q 20 -t species -r Bacteria @@write ${lib}_Bacteria.txt
#echo "Finished Bacterial taxonomy step."; endtax1=$(date +%T); echo "$endtax1"
mv ${lib}_Bacteria.txt post_taxonomy/

#get  taxonomic breakdown of kingdoms (l=1)
echo "Starting Domain step..."; starttax2=$(date +%T); echo "$starttax2"
paladin-plugins.py @@taxonomy -i $tsv -q 20 -t children -l 1 @@write ${lib}_Domains.txt
#echo "Finished Domain step."; endtax2=$(date +%T); echo "$endtax2"
mv ${lib}_Domains.txt post_taxonomy/

# get Eukaryota taxomomic breakdown
echo "Starting Bacterial taxonomy step..."; starttax1=$(date +%T); echo "$starttax1"
HTTP_PROXY=premise.sr.unh.edu:3128 paladin-plugins.py @@taxonomy -i $tsv -q 20 -t children -r Eukaryota @@write ${lib}_Eukaryota.txt
echo "Finished Bacterial taxonomy step."; endtax1=$(date +%T); echo "$endtax1"
mv ${lib}_Eukaryota.txt post_taxonomy/

#get  taxonomic breakdown of kingdoms (l=1)
echo "Starting phyla step..."; starttax2=$(date +%T); echo "$starttax2"
HTTP_PROXY=premise.sr.unh.edu:3128 HTTPS_PROXY=premise.sr.unh.edu:3128  paladin-plugins.py @@taxonomy -i $tsv -q 20 -t children -l 2 @@write ${lib}_phylaKingdom.txt
echo "Finished phyla step."; endtax2=$(date +%T); echo "$endtax2"
mv ${lib}_phylaKingdom.txt post_taxonomy/
fi
done
