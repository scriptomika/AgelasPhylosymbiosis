
#1 annotate full metagenome assembly with eggNOG-mapper
for assembly in  megahit_*/final.contigs.fa; do dir=${assembly%/*}; sp=${dir#*_};
do
#2 extract Metazon hits
cat ${dir}/${sp}.emapper.annotations | grep 'Metazoa' > ${dir}/${sp}.emapper.Metazoa
cut -f1 ${dir}/${sp}.emapper.Metazoa |awk -F'_' '{print $1"_"$2}' sort |uniq > contigs

#3 pull seqs from metagenome assembly
cat ${dir}/final.contigs.fa  |awk -F' ' '{print $1}' >test
selectSeqs.pl -f contigs test > ${dir}/sponge.final.contigs.fa

#4 reciprocalblast against Amphimedon
agelas="${dir}/sponge.final.contigs.fa"
aq="Amphimedon_queenslandica.faa"

blastx -query $agelas  -subject $aq -evalue 1e-50 -outfmt 6 -num_threads 24 -out blast1
tblastn -query $aq -subject $agelas -evalue 1e-50 -outfmt 6 -num_threads 24 -out blast2

./parse_recip_blast.py blast1 blast2 1  2 12 high rbh
tail -n+2 rbh |cut -f1 |sort|uniq >rbh2

#5 extract candidate sponge sequences
selectSeqs.pl -f rbh2 ${dir}/sponge.final.contigs.fa > ${dir}/sponge.final2.contigs.fa

#6 extract reads from sample fastqs to generate single-sample assemblies
brown="BZ076  KY039 KY116 KY141 KY142 KY177 SX011 SX012 SX031 SX057 SX112 KY174 KY196 KY197 KY198 BZ039 BZ055 BZ098 BZ099 BZ100 BZ305"
contub=
citrina=
sventres=
clathrodes=
mkdir sample_assemblies/

sample_assembly () {
   cut -f1 rbh2 | tr "\n" " " | xargs samtools view -bh ../../mapping2022_agelas/${1}.bam.sorted > ${1}_sponge.bam
   samtools sort -n ${1}_sponge.bam -o ${1}_sponge.bam.sorted
   samtools collate -u -O ${1}_sponge.bam.sorted | samtools fastq -1 ${1}_sponge1.fq -2 ${1}_sponge2.fq -0 /dev/null -s /dev/null -n 
   megahit -t 24 -o ${1}_megahit -1 ${1}_sponge1.fq -2 ${1}_sponge2.fq
   mv ${1}_megahit/final.contigs.fa sample_assemblies/${1}.final.contigs.fa
}
if [ $sp == "Acontub"]; then for lib in $contub; do sample_assembly $lib; done; fi
if [ $sp == "Abrown"]; then for lib in $brown; do sample_assembly $lib; done; fi
if [ $sp == "Acitrina"]; then for lib in $citrina; do sample_assembly $lib; done; fi
if [ $sp == "Asven"]; then for lib in $sventres; do sample_assembly $lib; done; fi
if [ $sp == "Aclath"]; then for lib in $clathrodes; do sample_assembly $lib; done; fi

done
