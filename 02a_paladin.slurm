#!/bin/bash

#SBATCH --partition=shared
#SBATCH -J AgPAL1
#SBATCH --output paladin1_1.log
#SBATCH --mem 500Gb
#SBATCH --open-mode=append
#SBATCH -x node117,node118
module load linuxbrew/colsa


# to run in parallel: distribute trimmed fastq across multiple dirs (trimmed_set1-6), duplicate this slurm 
# keep pairs separate for naming convenience. Run reads1 and reads2 separately

for fq in trimmed_set1/*1P.fastq
do 
libpath=${fq%.trim_*P.fastq}
lib=${libpath#trimmed_set1/}
echo $lib 			#msp39
echo $fq 			#msp39.trim_1P.fastq

#file1
paladin align /mnt/lustre/hcgs/shared/databases/uniref90/paladin/uniref90.fasta.gz \
$fq -t 24 -f 30 -o ${lib}_1_paladin -P http://premise.sr.unh.edu:3128

done

#file2, edit output string
for fq in trimmed_set1/*2P.fastq
do 
libpath=${fq%.trim_*P.fastq}
lib=${libpath#trimmed_set1/}
echo $lib 			#msp39
echo $fq 			#msp39.trim_2P.fastq

#file1
paladin align /mnt/lustre/hcgs/shared/databases/uniref90/paladin/uniref90.fasta.gz \
$fq -t 24 -f 30 -o ${lib}_2_paladin -P http://premise.sr.unh.edu:3128

done
#mv *.sam *tsv paladin_results

#alternative, running from login node if HTTP error:

# example:
# HTTP_PROXY=HTTPS_PROXY=premise.sr.unh.edu:3128 paladin-plugins.py @@uniprot -i msp39_2_paladin.sam -c "organism" "protein names" "genes" "pathway" "features" "go" "reviewed" "existence" "comments" "database(KEGG)" "database(GeneID)" "database(PATRIC)" "database(EnsemblBacteria)" "organism-id" "lineage(all)" -b 4000 @@write msp39_2_paladin_uniprot.tsv
