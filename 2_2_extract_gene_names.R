#!/usr/bin/Rscript

# extract PubMed genes that are related to PREGNANCY (based on abstract mining)

# 1) search the PubMed website with code words **tissue**  AND **gene**
# 2) download the text file
# 3) run Julius' script that that eliminates everything except abstracts
# 4) run "2_1..." script that further prunes abstracts
# 5) run this script which extracts valid gene names

#rm(list=ls())

# load the full collection human genes
hg=read.table("ucsc_HUGO_ENTREZ_chr1-23_withDescriptions_hg19_PROCESSED.txt",
              stringsAsFactors=F,h=F); dim(hg); head(hg); table(hg$V1) # ~/Biostuff/hg19_HUMAN_GENES/
colnames(hg)=c("CHR","START","END","ENTREZ","HUGO","Description")
hg_genes=hg[which(nchar(hg$HUGO)>1),"HUGO"]; length(hg_genes) # 2 = arguably too much risk with short gene names
#hg[grep("IL+[1-2]{1,1}$",hg$HUGO),][1:10,]

# explore how many letters usually are there in the HUGO gene name (from left side)
#ltrs = NULL; for (i in 1:dim(hg)[1]) ltrs = c(ltrs, unlist(strsplit(hg[i,"HUGO"],"[0-9]"))[1])
#head(ltrs); length(ltrs); dim(hg)
##hist(nchar(ltrs),col="grey")
#ltrs[which(nchar(ltrs)>8)]
#hg[which(nchar(ltrs)==1),1:5]

# load dangerous gene names that are also biomed acronyms
#acronym_dir="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY/"
acronyms_tbl = read.table("restricted_acronyms-geneNames.txt",sep="\t",h=T,stringsAsFactors=F)
restricted_acronyms = sort(unique(acronyms_tbl$Acronym))
# majority of acronyms are present in the gene-name TRANSLATOR file..
# .. and thus can be detected via their "long-name"
#hg[which(hg$HUGO %in% restricted_acronyms),]


# load what was generated in previous script (cleaned abstracts)
out_dir="./PubMed_PRUNE/" #~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY
load(paste(out_dir,"cleaned_abstracts_GYNECOLOGY.RData",sep=""))


phenotypes = names(cleaned_abstracts_exclusivityON)  # assumption!  both files have the same phenotypes
phenotypes = phenotypes[-grep("stats|HEART",phenotypes)]

gene_lists = list() # here tables of gene-freqs for all phenotypes will be accumulated

# cycle through various combinations of phenotypes/methods
for (phenotype in phenotypes) {
        for (exclusivity in c("xcl","nxc")) {
                for (translation in c("trn","unt")) {
                        
                        # report
                        obj_name = paste(phenotype,exclusivity,translation,sep="_")
                        print(obj_name)
                        
        # get the relevant data
        if (exclusivity=="xcl") tempor = get("cleaned_abstracts_exclusivityON")
        if (exclusivity=="nxc") tempor = get("cleaned_abstracts_exclusivityOFF")
        raw.txt = tempor[[phenotype]]; length(raw.txt); rm(tempor)
        
        # DECISION WHETHER TRANSLATOR should be used
        if (translation=="trn") {
                #transl_dir="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY/"
                translator = read.table("TRANSLATOR_misspelled_gene_names.txt",
                                        stringsAsFactors=F,h=T,sep="\t")# TRANSLATOR OF SOME (!) GENE NAMES
                from_length = nchar(translator$from)
                translator = translator[order(from_length,decreasing = T),]
        } else { translator = data.frame(from="111111",to="222222",stringsAsFactors = F) }

        
         
cumm1=NULL  # cummulation of potential gene names extracted from abstracts without TRANSLATOR
cumm2=NULL  # cummulation of REAL gene names extracted from abstracts using TRANSLATOR
cumm3=NULL  # cummulation of ALL unique-in-one-abstract gene names (ddtected with or without translator). for "times-mentioned" threshold

### perform cycling through abstracts with gene-name extraction procedure
for (j in 1:length(raw.txt)) { # takes a while... quite a big one
                
        #print(paste(j,"/",length(raw.txt),sep=""))
        # this is the text that we will be working with in this cycle
        txt=raw.txt[j]
        
        ##  use the TRANSLATOR to EXTARCT (!) gene names that are written in a non-standard manner ("long-format")
        replaced = NULL
        for (t in 1:dim(translator)[1]) { # note that the order of translator is important. longer entries first!
                regexp_phrase = paste("([[:punct:]]|\\s)+", translator[t,"from"] ,"([[:punct:]]|\\s)+",sep="")
                # note that currently IL-XX genes are translated problematically (in case of receptr - is assignes erroneously)
                if (length(grep(regexp_phrase,txt))>0) {
                        replaced = c(replaced, translator[t,"to"])
                        txt = paste(unlist(strsplit(txt, translator[t,"from"])),collapse=" REPLACED ")
                }
                rm(regexp_phrase)
        }

        # remove all possible punctuation marks, but preserve the hyphen
        txt=paste(unlist(strsplit(txt,"-")),collapse="zzzzz") # "zubiquitilation of hyphens"
        txt=paste(unlist(strsplit(txt,"[[:punct:]]")),collapse=" ")
        txt=paste(unlist(strsplit(txt,"zzzzz")),collapse="-")
        
        # extract all remaining POTENTIAL gene names (assumptions are made below!)
        lst1=unlist(regmatches(txt, gregexpr('([[:punct:]]|[[:upper:]]|[0-9])+', txt)))
        lst2=unlist(regmatches(lst1, gregexpr('^[[:upper:]]{2,20}.*', lst1)))
        lst2 = sort(unique(lst2))
        lst2 = lst2[which(! lst2 %in% "REPLACED")]  # do not include artifacts

        # ELIMINATE GENES THAT ARE EASILY CONFUSED WITH THE TERMS FROM ..
        # .. OBSTETRIC/MEDICAT/TECHNICAL FIELDS. eliminate only from the first set!
        lst3 = lst2[which(! lst2 %in% restricted_acronyms)]
        
        cumm1=c(cumm1,unique(lst3))
        cumm2=c(cumm2,unique(replaced))
        cumm3=c(cumm3,unique(c(lst3,replaced))) # will be used to estimate number of instances when a gene was mentioned
        rm(txt,lst1,lst2,lst3,replaced)
}


# directly extracted POTENTIAL genes
genes_1 = cumm1[which(cumm1 %in% hg_genes)]
#hist(as.numeric(table(genes_1)),breaks=100,col="grey")
#genes_1 = sort(unique(genes_1))

# genes detected/extracted via TRANSLATOR file
genes_2 = cumm2[which(cumm2 %in% hg_genes)]
#hist(as.numeric(table(genes_2)),breaks=100,col="grey")
#genes_2 = sort(unique(genes_2))

# detected via BOTH METHODS
genes_3 = cumm3[which(cumm3 %in% hg_genes)]
#hist(as.numeric(table(genes_3)),breaks=100,col="grey")
#genes_3 = sort(unique(genes_3))

#library(gplots)
#temp = list(extracted = genes_1,translated = genes_2)
#venn(temp); rm(temp)

# all together
#if (length(genes_2)==0) {genes = sort(genes_1)} else {genes=sort(c(genes_1,genes_2))}
genes = genes_3

temp1=data.frame(gene=names(table(genes)),freq=as.numeric(table(genes)),stringsAsFactors = F)
temp2= merge(temp1,hg,by.x="gene",by.y="HUGO",all.x=T); rm(temp1)
gene.freq=temp2[rev(order(temp2$freq)),]; rm(temp2)
#gene.freq[1:20,]

obj_name = paste(phenotype,exclusivity,translation,sep="_")
gene_lists[[obj_name]]=gene.freq
rm(gene.freq, genes, genes_1,genes_2, cumm1,cumm2,cumm3,obj_name)

} # end of cycling through various phenotypes
}
}


# also store the version (version-control) identificator (hash)
hash = system("git log --pretty=format:'%h' -n 1",intern=TRUE)
gene_lists[["hash"]] = hash
#gene_lists$hash

# save the R object with results 
result_dir = "./PubMed_GENES/" #~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY/
save(list=c("gene_lists"),file=paste(result_dir,"PubMed_extracted_genes_GYNECOLOGYandCONTROL.RData",sep=""))

quit()


############################################
###########   PREVIEW

library(gplots)

# for pregnancy-related genes
temp=list(endom.=gene_lists$ENDOMETRIUM$gene,myom.=gene_lists$MYOMETRIUM$gene,
          cervix = gene_lists$CERVIX$gene,uterus = gene_lists$UTERUS$gene,
          placenta = gene_lists$PLACENTA$gene)
venn(temp)

# for control-set of tissues/phenotypes
temp=list(bladder = gene_lists$BLADDER$gene,bone=gene_lists$BONE$gene,
          penile = gene_lists$PENILE$gene,prostate = gene_lists$PROSTATE$gene,
          trachea = gene_lists$TRACHEA$gene)
venn(temp)

# for control-set of tissues/phenotypes

temp=list(intest = gene_lists$INTESTINE_xcl_trn$gene,liver=gene_lists$LIVER_xcl_trn$gene,
          lungs = gene_lists$LUNGS_xcl_trn$gene,skin = gene_lists$SKIN_xcl_trn$gene,
           muscle= gene_lists$MUSCLE_xcl_trn$gene)
venn(temp)







