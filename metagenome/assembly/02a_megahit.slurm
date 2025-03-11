#!/bin/bash
#SBATCH --partition=shared
#SBATCH -J AassAbrown
#SBATCH --output 02A_megahit.Abrown.log
#SBATCH --open-mode=append
#SBATCH --mem=512000

module load linuxbrew/colsa

# SYNOPSIS
# Run MegaHit to assembly metagenome contigs for each Agelas genetic clade
# (conifera/tubulata; dispar/brown; clathrodes, sventres, citrina)
# break up to run in parallel

#dispar/brown
for lib in BZ076 BZ305 KY039 KY116 KY141 KY142 KY177 SX011 SX012 SX031 SX057 SX112 KY174 KY196 KY197 KY198
do 
echo $lib
cat ../trimmed_done/${lib}*1P.fastq >> Abrown.trim_1P.fastq
cat ../trimmed_done/${lib}*2P.fastq >> Abrown.trim_2P.fastq
done

megahit -t 24 -o megahit_Abrown \
-1 Abrown.trim_1P.fastq -2 Abrown.trim_2P.fastq 


#citrina
for lib in CU22 CU24 CU26 CU30 KY122 KY175 KY186 KY202 KY209 SX134 SX188 SX190
do 
echo $lib
cat ../trimmed_done/${lib}*1P.fastq >> Acitrina.trim_1P.fastq
cat ../trimmed_done/${lib}*2P.fastq >> Acitrina.trim_2P.fastq
done

megahit -t 24 -o megahit_Acitrina \
-1 Acitrina.trim_1P.fastq -2 Acitrina.trim_2P.fastq 


#clathrodes
for lib in BZ142 KY120 KY121 KY123 KY144 KY201 SX033 SX135 SX149 SX191 SX311
do 
echo $lib
cat ../trimmed_done/${lib}*1P.fastq >> Aclath.trim_1P.fastq
cat ../trimmed_done/${lib}*2P.fastq >> Aclath.trim_2P.fastq
done

megahit -t 24 -o megahit_Aclath \
-1 Aclath.trim_1P.fastq -2 Aclath.trim_2P.fastq 

#sventres
for lib in BZ075 BZ231 BZ290 BZ291 BZ303 CU002 CU003 CU004 CU005 CU008 SX002 SX003 SX058 SX192 SX309
do 
echo $lib
cat ../trimmed_done/${lib}*1P.fastq >> Asven.trim_1P.fastq
cat ../trimmed_done/${lib}*2P.fastq >> Asven.trim_2P.fastq
done

megahit -t 24 -o megahit_Asven \
-1 Asven.trim_1P.fastq -2 Asven.trim_2P.fastq 

#conifera/tubulata
for lib in BZAgtub17 BZAgtub20 CUAgc15 CUAgc16 msp25 msp26 msp27 msp28 msp30 msp31 msp32 msp35 msp43 msp44 msp45 SXAgco05 SXAgco13 KYAgtu01 KY100Agt1 KY100Agt2 KY100Agt3 KY100Agt4 KY150Agt1 KY150Agt3 KY150Agt2 KY150Agt4 KY150Agt5 KY200Agt1 KY200Agt2 KY200Agt3 KY200Agt4 KY200Agt5
do 
echo $lib
cat ../trimmed_done/${lib}*1P.fastq >> Acontub.trim_1P.fastq
cat ../trimmed_done/${lib}*2P.fastq >> Acontub.trim_2P.fastq
done

megahit -t 24 -o megahit_Acontub \
-1 Acontub.trim_1P.fastq -2 Acontub.trim_2P.fastq 
