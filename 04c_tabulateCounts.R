#!/bin/Rscript
#usage Rscript --vanilla ./04c_tabulate_counts.R
#takes all files ending in "_counts.txt" (produced from pathway-specific paladin UniProt search)
# tabulates to produce long table (CSV) with 
# library, taxon domain (metazoa, bacteria, archea), ec code, path name, total counts mapped

#set pathway to count files
setwd("./post_pathways/counts/")

loopthrufiles<-function(countfile)
{
print(countfile)
lib<-unlist(strsplit(countfile,split="_"))[1]
ec<-unlist(strsplit(countfile,split="_"))[3]
domain<-unlist(strsplit(countfile,split="_"))[2]
df<-read.table(countfile,row.names=1,skip=1) #skip header line since contains spaces
df<-df[1:dim(df)[1]-1,]
tot<-rowSums(df)
df.out<-data.frame(lib=rep(lib,length(tot)),domain=rep(domain,length(tot)),ec=rep(ec,length(tot)),enz=rownames(df),totalcount=tot)
return(df.out)
}

#files<-list.files(pattern="_counts.txt")
files<-intersect(list.files(pattern="BZ"),list.files(pattern="counts.txt"))
dflist<-lapply(files, loopthrufiles)
df.save <- do.call("rbind", dflist)
rownames(df.save)<-seq(1:nrow(df.save))
write.csv(df.save, "pathway_enzyme_counts.BZ.csv")

files<-intersect(list.files(pattern="KY"),list.files(pattern="counts.txt"))
dflist<-lapply(files, loopthrufiles)
df.save <- do.call("rbind", dflist)
rownames(df.save)<-seq(1:nrow(df.save))
write.csv(df.save, "pathway_enzyme_counts.KY.csv")

files<-intersect(list.files(pattern="SX"),list.files(pattern="counts.txt"))
dflist<-lapply(files, loopthrufiles)
df.save <- do.call("rbind", dflist)
rownames(df.save)<-seq(1:nrow(df.save))
write.csv(df.save, "pathway_enzyme_counts.SX.csv")

files<-intersect(list.files(pattern="CU"),list.files(pattern="counts.txt"))
dflist<-lapply(files, loopthrufiles)
df.save <- do.call("rbind", dflist)
rownames(df.save)<-seq(1:nrow(df.save))
write.csv(df.save, "pathway_enzyme_counts.CU.csv")

files<-intersect(list.files(pattern="msp"),list.files(pattern="counts.txt"))
dflist<-lapply(files, loopthrufiles)
df.save <- do.call("rbind", dflist)
rownames(df.save)<-seq(1:nrow(df.save))
write.csv(df.save, "pathway_enzyme_counts.msp.csv")
