#!/bin/bash

#SBATCH --partition=shared
#SBATCH -J Trim
#SBATCH --output trim.log
#SBATCH --open-mode=append

module load linuxbrew/colsa

for fq in *_1.fq.gz

do 
echo $fq
lib=${fq%%_*};
echo $lib;

trimmomatic PE -threads 24 \
-baseout ${lib}.trim.fastq \
${lib}_*_1.fq.gz \
${lib}_*_2.fq.gz \
LEADING:3 TRAILING:3 \
ILLUMINACLIP:nebnext.fa:2:30:10 MINLEN:55

done
