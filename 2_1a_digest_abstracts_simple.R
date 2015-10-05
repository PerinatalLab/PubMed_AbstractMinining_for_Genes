#!/usr/bin/Rscript

## This is an alternative script to Jonas's 2_1_digest_abstracts.R.
## This version is simplified and includes no special pruning. It is also assumed
## that all required phenotypes have been combined into one file -
## i.e., in this script files are processed and exported one-by-one.

# 1) search the PubMed website with code words **tissue**  AND **gene**
# 2) download the text file
# 3) run Julius' script that that eliminates everything except abstracts
# 4) run this script for some basic cleanup

## USAGE: ./2_1a_digest_abstracts_simple.R path_to_dir_with_PubMed_subdirs

args=commandArgs(TRUE)
working_dir=args[1]

PubMedDir=paste(working_dir,"PubMed_DIGEST/",sep="")
out_dir=paste(working_dir,"PubMed_PRUNE/",sep="")
system(paste("mkdir ",out_dir,sep=""))

file_list = list.files(PubMedDir)
files_ok = file_list[grep("abstracts",file_list)]

for(file_name in files_ok){
  print(paste("working on file",file_name))
  raw.txt=readLines(paste(PubMedDir,file_name,sep=""))
  print(paste("number of abstracts (initial): ",length(raw.txt),sep=""))
  n_abstracts_0 = length(raw.txt)

  goo_length = which(nchar(raw.txt)>100)
  raw.txt=raw.txt[goo_length]; rm(goo_length)
  print(paste("number of abstracts (length > 100 smbls): ",length(raw.txt),sep=""))
  n_abstracts_1 = length(raw.txt)

  # get rid of tab symbol in the begining of the text string
  for (i in 1:length(raw.txt)) raw.txt[i] = unlist(strsplit(raw.txt[i],"\t"))[2]

  name_chunks = unlist(strsplit(file_name,"_"))
  name_chunks[3] = "PRUNED"
  new_name = paste(name_chunks,sep="_",collapse="_")
  write.table(raw.txt,paste(out_dir,new_name,sep=""),quote = F,row.names = F,col.names = F)
}