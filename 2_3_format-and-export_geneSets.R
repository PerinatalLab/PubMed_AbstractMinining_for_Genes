#!/usr/bin/Rscript

########################################################################
########################################################################
######################### export the gene sets

rm(list=ls())  # cleanup

### option for the DELIVERY-related of phenotypes ***
result_dir = "./PubMed_GENES/" #~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jul_OBSTETRICS
load(paste(result_dir,"PubMed_extracted_genes_OBSTETRICSandCONTROL.RData",sep=""))

obj_lst = names(gene_lists)
obj_lst = obj_lst[-grep("hash",obj_lst)]

inrich_input = NULL
for (z in 1:length(obj_lst)) {  # for each type of set/settings
        gene.freq = gene_lists[[obj_lst[z]]]
        col1 = data.frame(ENTREZ = gene.freq[gene.freq$freq>=1,"ENTREZ"],
                                   GENESET=paste("PM",":",obj_lst[z],"_1",sep=""),
                                   Descript=".")
        col2 = data.frame(ENTREZ = gene.freq[gene.freq$freq>=2,"ENTREZ"],
                                   GENESET=paste("PM",":",obj_lst[z],"_2",sep=""),
                                   Descript=".")
        col=rbind(col1,col2)
        inrich_input = rbind (inrich_input,col)
        rm(col,col1,col2,gene.freq)
}

#dim(inrich_input)
#table(inrich_input$GENESET)
#head(inrich_input)
 
colnames(inrich_input)[1]=paste("##",colnames(inrich_input)[1],sep="")
folder = "./PubMed_GENES/" #~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jul_OBSTETRICS
file_name = paste(folder,"PubMed_geneSets_forINRINCH_OBSTETRICSandCONTROL.txt",sep="")
write.table(inrich_input,file_name,row.names=F,col.names=T,sep="\t",quote=F)
        