
########################################################################
########################################################################
######################### export the gene sets

rm(list=ls())  # cleanup

### option for the pregnancy-related of phenotypes(tissues) ***
load("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/obgn_genes.RData")
load("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/cntrl_genes.RData")
lst=ls()
obj_lst = lst[grep("^obg_.{7,7}$|^ctrl_.{7,7}$",lst)]  # gene frequencies per each phenotype/tissue

### option for the control-set of phenotypes(tissues) ***
#load("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/cntrl_genes.RData")
#lst=ls()
#obj_lst = lst[grep("^ctrl_.{7,7}$|^ctrl_.{7,7}$",lst)]  # gene frequencies per each phenotype/tissue

COLLECTOR = NULL
for (z in 1:length(obj_lst)) {  # for each type of settings
        temp_obj = get(obj_lst[z])
        types = names(temp_obj)
        
        # optional (depending on whether we want to include everything or only relevant)   ***
        types = types[ which( !types %in% c("PLACENTA","DENTAL","TRACHEA"))]  # suppress or leave 
        
        collector = NULL
        for (type in types) { # for each type of tissue
                gene.freq=temp_obj[[type]]        
                temp1 = data.frame(ENTREZ = gene.freq[gene.freq$freq>=1,"ENTREZ"],
                                   GENESET=paste("PM",":",type,"_",paste(unlist(strsplit(obj_lst[z],
                                        "_"))[2:3],collapse="_"),"_1",sep=""),Descript = ".")
                
                temp2 = data.frame(ENTREZ = gene.freq[gene.freq$freq>=2,"ENTREZ"],
                                   GENESET=paste("PM",":",type,"_",paste(unlist(strsplit(obj_lst[z],
                                        "_"))[2:3],collapse="_"),"_2",sep=""),Descript = ".")
                collector=rbind(collector,temp1,temp2); rm(gene.freq,temp1,temp2)
        }
        COLLECTOR = rbind(COLLECTOR,collector)
        rm(collector, types, temp_obj)
}

        folder = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_GENES/"
        colnames(COLLECTOR)[1]=paste("##",colnames(COLLECTOR)[1],sep="")
        file_name = paste(folder,"PubMed_geneSets_forINRICH_GYNECOLOGYandCONTROL.txt",sep="")
        write.table(COLLECTOR,file_name,row.names=F,col.names=T,sep="\t",quote=F)
        rm(file_name,folder)


