
########################################################################
########################################################################
######################### export the gene sets

rm(list=ls())  # cleanup

### option for the DELIVERY-related of phenotypes ***
result_dir = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/" 
load(paste(result_dir,"PubMed_extracted_genes_DELIVERY.RData",sep=""))

prtr.trn = gene_lists$PARTURITION_trn
prtr.unt = gene_lists$PARTURITION_unt
obj_lst = c("prtr.trn","prtr.unt")

for (z in 1:length(obj_lst)) {  # for each type of settings
        gene.freq=get(obj_lst[z])
        col1 = data.frame(ENTREZ = gene.freq[gene.freq$freq>=1,"ENTREZ"],
                                   GENESET=paste("PM",":",obj_lst[z],sep=""),
                                   Descript="obs>=1")
                
        col2 = data.frame(ENTREZ = gene.freq[gene.freq$freq>=2,"ENTREZ"],
                                   GENESET=paste("PM",":",obj_lst[z],sep=""),
                                   Descript="obs>=2")
        rm(gene.freq)
        
        folder = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_GENES/"
        
        
        # no restriction on gene frequency
        printout=col1
        colnames(printout)[1]=paste("##",colnames(printout)[1],sep="")
        file_name = paste(folder,"PubMed_",obj_lst[z],".min1","_xxxxxx.txt",sep="")
        write.table(printout,file_name,row.names=F,col.names=T,sep="\t",quote=F)
        
        
        printout=col2
        colnames(printout)[1]=paste("##",colnames(printout)[1],sep="")
        file_name = paste(folder,"PubMed_",obj_lst[z],".min2","_xxxxxx.txt",sep="")
        write.table(printout,file_name,row.names=F,col.names=T,sep="\t",quote=F)
        rm(col1,col2,file_name,printout,folder)
        
        }
        




