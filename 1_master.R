#!/usr/bin/Rscript

# passes arguments to the bash script that takes care of PubMed abstracts

setwd("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY")
bash_script = "./1_extract_Abstracts_from_PubMed_output.sh"

# abstracts that are downloaded from PubMed
file_list1 = list.files("./PubMed_RAW/")
file_list = file_list1[grep("SKIN|INTEST|MUSCL|HEART|LIVER|LUNG",file_list1)]


# naming convention:  PubMed_webSearch_RAW_CERVIX_GENES_2014Jun17_n1362.txt
for (i in 1:length(file_list)) {
file_name = file_list[i]
print(file_name)
name_chunks = unlist(strsplit(file_name,"_"))
name_chunks[3] = "DIGESTED"
name_chunks[7] = unlist(strsplit(name_chunks[7],"\\."))[1] # get rid of .txt
new_name = paste(name_chunks,sep="_",collapse="_")
cmnd = paste(bash_script,file_name,new_name,sep=" ")
system(cmnd,ignore.stdout = F)
}



