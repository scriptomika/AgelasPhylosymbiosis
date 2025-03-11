#!/bin/Rscript
#usage Rscript --vanilla ./04d_pathway_contributions.R
#takes all files ending in "_pathways.csv"

#produced from paladin @@pathway output modified:
# for f in *pathways.txt; do cut -f2 --complement $f |perl -p -e 's/\t/,/g' > ${f%*.txt}.csv

# tabulates to produce long table (CSV) with 
# library, ec code, taxon, proportion of pathway enzymes recovered in reads from taxon

library(tidyr)
#set pathway to pathway output files
setwd("./post_pathways/")

loopthrufiles<-function(contribfile)
{
print(contribfile)
lib<-unlist(strsplit(contribfile,split="_"))[1]
df<-read.csv(contribfile,row.names=1,header=T)  %>% tibble::rownames_to_column(var="ec") %>% tidyr::pivot_longer(.,cols=starts_with("X1"),names_to="taxon",values_to="contrib")
df$taxon<-gsub("X1.","",df$taxon)
df$lib<-rep(lib,nrow(df))
df.out<-data.frame(df)
return(df.out)
}

files<-list.files(pattern="_pathways.csv")
dflist<-lapply(files, loopthrufiles)

df.save <- do.call("rbind", dflist)
rownames(df.save)<-seq(1:nrow(df.save))
write.csv(df.save, "pathway_contributions_tabulated.txt")
