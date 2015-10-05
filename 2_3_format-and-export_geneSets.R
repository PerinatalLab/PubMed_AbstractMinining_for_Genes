#!/usr/bin/Rscript

## This script reformats and exports the gene sets, also filtering by the min number of mentions

## USAGE: ./2_3_format-and-export_geneSets.R working_dir/, minimum_mentions

options(stringsAsFactors = F)
args=commandArgs(TRUE)
working_dir=args[1]
min_mentions=args[2] ## genes with lower number of mentions will not be included

gene_dir = paste(working_dir,"PubMed_GENES/",sep="")
file_list = list.files(gene_dir)
file_list = file_list[grep(".txt$",file_list)]

########## simple version

for(file_name in file_list){
    print(paste("working on file",file_name))
    gene.freq=read.table(paste(gene_dir,file_name,sep=""),h=F)
    out_name = paste(gene_dir,unlist(strsplit(file_name,"\\."))[1],".set",sep="")
    temp1 = data.frame(ENTREZ = gene.freq[gene.freq$V2>=min_mentions,6], GENESET="NASAL_POLYP")
    write.table(temp1,out_name,quote = F,row.names = F,col.names = F)
}
