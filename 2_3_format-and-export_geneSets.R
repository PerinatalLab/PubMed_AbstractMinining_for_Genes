
########################################################################
########################################################################
######################### export the gene sets

rm(list=ls())  # cleanup

### option for the pregnancy-related of phenotypes(tissues) ***
load("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/obgn_genes.RData")
lst=ls()
obj_lst = lst[grep("^obg_.{7,7}$|^ctrl_.{7,7}$",lst)]  # gene frequencies per each phenotype/tissue

### option for the control-set of phenotypes(tissues) ***
#load("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/cntrl_genes.RData")
#lst=ls()
#obj_lst = lst[grep("^ctrl_.{7,7}$|^ctrl_.{7,7}$",lst)]  # gene frequencies per each phenotype/tissue

for (z in 1:length(obj_lst)) {  # for each type of settings
        temp_obj = get(obj_lst[z])
        types = names(temp_obj)
        
        # optional (depending on whether we want to include everything or only relevant)   ***
        types = types[ which( !types %in% c("PLACENTA","DENTAL","TRACHEA"))]  # suppress or leave 
        
        col1 = col2 = NULL  # collectors
        for (type in types) { # for each type of tissue
                gene.freq=temp_obj[[type]]        
                temp1 = data.frame(ENTREZ = gene.freq[gene.freq$freq>=1,"ENTREZ"],
                                   GENESET=paste("PM",":",type,sep=""),
                                   Descript=paste("obs>=1",paste(unlist(strsplit(obj_lst[z],
                                                                                 "_"))[2:3],collapse=","),sep=","))
                col1=rbind(col1,temp1)
                temp2 = data.frame(ENTREZ = gene.freq[gene.freq$freq>=2,"ENTREZ"],
                                   GENESET=paste("PM",":",type,sep=""),
                                   Descript=paste("obs>=2",paste(unlist(strsplit(obj_lst[z],
                                                                                 "_"))[2:3],collapse=","),sep=","))
                col2=rbind(col2,temp2); rm(gene.freq,temp1,temp2)
        }
        
        middle = paste(unlist(strsplit(obj_lst[z],"_")),collapse=".")
        folder = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_GENES/"
        
        # no restriction on gene frequency
        printout=col1
        colnames(printout)[1]=paste("##",colnames(printout)[1],sep="")
        file_name = paste(folder,"PubMed_",middle,".min1","_",length(types),"pheOvrlp.txt",sep="")
        write.table(printout,file_name,row.names=F,col.names=T,sep="\t",quote=F)
        
        printout=col2
        colnames(printout)[1]=paste("##",colnames(printout)[1],sep="")
        file_name = paste(folder,"PubMed_",middle,".min2","_",length(types),"pheOvrlp.txt",sep="")
        write.table(printout,file_name,row.names=F,col.names=T,sep="\t",quote=F)
        rm(col1,col2,file_name,printout,middle,folder,temp_obj,types)
}


### below need a review




