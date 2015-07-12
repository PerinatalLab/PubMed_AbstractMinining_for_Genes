
########################################################################
########################################################################
######################### export the gene sets

rm(list=ls())  # cleanup

### option for the pregnancy-related of phenotypes(tissues) ***
load("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY/WORK_FILES/PubMed_extracted_genes.RData")

names(gene_lists)

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

        folder = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY/PubMed_GENES/"
        colnames(collector)[1]=paste("##",colnames(collector)[1],sep="")
        file_name = paste(folder,"PubMed_geneSets_forINRICH_GYNECOLOGYandCONTROL.txt",sep="")
        write.table(collector,file_name,row.names=F,col.names=T,sep="\t",quote=F)
        rm(file_name,folder)


