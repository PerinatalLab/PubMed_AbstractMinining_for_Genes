#!/usr/bin/Rscript

########################################################################
########################################################################
######################### export the gene sets

rm(list=ls())  # cleanup

load("./PubMed_GENES/PubMed_extracted_genes.RData") #~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY

types = names(gene_lists)

collector = NULL
for (type in types) {  # for each type of settings

        gene.freq = gene_lists[[type]]
        temp1 = data.frame(ENTREZ = gene.freq[gene.freq$freq>=1,"ENTREZ"],
                                   GENESET=paste("PM",":",type,"_1",sep=""),Descript = ".")
        temp2 = data.frame(ENTREZ = gene.freq[gene.freq$freq>=2,"ENTREZ"],
                                   GENESET=paste("PM",":",type,"_2",sep=""),Descript = ".")
                collector=rbind(collector,temp1,temp2); rm(gene.freq,temp1,temp2)
}
dim(collector)
head(collector)

        folder = "./PubMed_GENES/" #~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY
        colnames(collector)[1]=paste("##",colnames(collector)[1],sep="")
        file_name = paste(folder,"PubMed_geneSets_forINRICH_GYNECOLOGYandCONTROL.txt",sep="")
        write.table(collector,file_name,row.names=F,col.names=T,sep="\t",quote=F)
        rm(file_name,folder)


