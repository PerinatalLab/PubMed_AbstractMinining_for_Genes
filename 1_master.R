#!/usr/bin/Rscript

# passes arguments to the bash script that takes care of PubMed abstracts


bash_script = "./1_extract_Abstracts_from_PubMed_output.sh"

# abstracts that are downloaded from PubMed
file_list = list.files("./PubMed_RAW/")

for (i in 1:length(file_list)) {
        file_name = file_list[i]
        print(file_name)
        name_chunks = unlist(strsplit(file_name,"[[:punct:]]"))[1:2]
        new_name = paste(paste(name_chunks,sep="_",collapse="_"),"DIGESTED",sep="_")
        print(new_name)
        cmnd = paste(bash_script,file_name,new_name,sep=" ")
        system(cmnd,ignore.stdout = F)
}



