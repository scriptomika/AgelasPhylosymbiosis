#!/bin/bash

#run from shell on login due to HTTP errors on cluster nodes
for tsv in post_aggreg/*.both_uniprot.tsv

do 

libpath=${tsv%%.both_uniprot.tsv}
lib=${libpath#post_aggreg/}
echo $lib
if [ ! -f post_pathways/${lib}* ]
then
# report metabolic pathway 
echo "Starting pathways step..."; startPath=$(date  +%T); echo "$startPath"
#HTTPS_PROXY=HTTP_PROXY=premise.sr.unh.edu:3128 paladin-plugins.py @@pathways -i $tsv -q 20 -l 2 @@write ${lib}_pathways.txt
#try export proxy above instead
paladin-plugins.py @@pathways -i $tsv -q 20 -l 2 @@write ${lib}_pathways.txt


echo "Finished pathways step."; endPath=$(date  +%T); echo "$endPath"
mv ${lib}_pathways.txt post_pathways/
fi
done
