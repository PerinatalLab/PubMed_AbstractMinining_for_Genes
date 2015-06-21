
# extract PubMed genes that are related to PREGNANCY (based on abstract mining)

# 1) search the PubMed website with code words **tissue**  AND **gene**, download the text file
# 2) run Julius' script that that eliminates everything except abstracts
# 3) run this script that further prunes abstracts and extracts gene names

#rm(list=ls())

# load the full collection human genes
hg=read.table("~/Biostuff/hg19_HUMAN_GENES/ucsc_HUGO_ENTREZ_chr1-23_withDescriptions_hg19_PROCESSED.txt",
              stringsAsFactors=F,h=F); dim(hg); head(hg); table(hg$V1)
colnames(hg)=c("CHR","START","END","ENTREZ","HUGO","Description")
hg_genes=hg[which(nchar(hg$HUGO)>1),"HUGO"]; length(hg_genes) # 2 = arguably too much risk with short gene names
#hg[grep("IL+[1-2]{1,1}$",hg$HUGO),][1:10,]


# explore how many letters usually are there in the HUGO gene name (from left side)
#ltrs = NULL; for (i in 1:dim(hg)[1]) ltrs = c(ltrs, unlist(strsplit(hg[i,"HUGO"],"[0-9]"))[1])
#head(ltrs); length(ltrs); dim(hg)
##hist(nchar(ltrs),col="grey")
#ltrs[which(nchar(ltrs)>8)]
#hg[which(nchar(ltrs)==1),1:5]


###  get typical disease-names from definitions in ICD-code-file  (in order to get rid of abstracts with other diseases)
txt = read.table("~/Biostuff/ICD-9-CM-v32-master-descriptions/CMS32_DESC_SHORT_DX.txt",h=F,sep="\t",stringsAsFactors = F)
all_words=NULL  # takes about 5 seconds
for (i in 1:dim(txt)[1]) {
        lst1=unlist(regmatches(txt[i,1], gregexpr("[A-Za-z]+", txt[i,1])))
        all_words=c(all_words,lst1[which(nchar(lst1)>=3)])
        rm(lst1)
}
unq_words = sort(unique(all_words))
diseases1 = unq_words[grep("itis$|emia$|enia$|osis$|asia$|oma$|esis$|oses$|asis$|ydia$|rrhea$|onia$",unq_words)]
diseases2 = tolower(diseases1) # "ignore-case" option in grep will be used anyway
diseases = sort(unique(diseases2))
length(diseases)
#table(nchar(diseases))
#diseases[which(nchar(diseases)<6)]

##  get the list of animal-related terms that should not be in the abstracts (non-human subjects)
anim1 = read.table("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/animal_name_indicators.txt",
                   stringsAsFactors = F,h=F,sep="\t")
anim1=sort(unique(tolower(anim1[,1])))
anim2 = unique(c("animal","cattle","buffalo","ruminant", "cow","dog","rat","pig","cat","lion","horse","monkey","mouse",
                 "cows","dogs","rats","pigs","cats","lions","horses","monkeys","mice"))
animals=sort(unique(c(anim1,anim2))); rm(anim1,anim2)


PubMedDir="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_DIGEST/"
#file_list = list.files(PubMedDir,pattern="PLACEN|CERVIX|MYOMETR|ENDOMETR|UTER")  # pregnancy-related genes  ***
file_list = list.files(PubMedDir,pattern="PENILE|BLADD|BONE|DENTAL|PROSTAT|TRACHE")  # control set of genes (other tissues) ***
files_ok = file_list[grep("abstracts",file_list)]
pheno=NULL; for (i in 1:length(files_ok))pheno = c(pheno, unlist(strsplit(files_ok[i],"_"))[4]); print(pheno)

#phes = c("PREGNANCY","CERVIX","ENDOMETRIUM","MYOMETRIUM","PRETERM","BIRTH","GESTATIONAL","UTERINE","PLACENTA")
#phes =  c("ENDOMETRIUM","MYOMETRIUM","CERVIX","UTERUS","PLACENTA")
phes = pheno
        
gene_lists = list() # here tables for all phenotypes will be accumulated
summary_lists = list() # here will be collected numbers of abstracts remaining after each step

exclusivity_pruning = FALSE
translator_usage = TRUE

for (phe in phes) {
        print(phe)
        file_name = files_ok[which(pheno==phe)]
        raw.txt=readLines(paste(PubMedDir,file_name,sep=""))
        print(paste("number of abstacts (initial): ",length(raw.txt),sep=""))
        n_abstracts_0 = length(raw.txt)
        
        goo_length = which(nchar(raw.txt)>100)
        raw.txt=raw.txt[goo_length]; rm(goo_length)
        print(paste("number of abstacts (length > 100 smbls): ",length(raw.txt),sep=""))
        n_abstracts_1 = length(raw.txt)
        
        # get rid of tab symbol in the begining of the text string
        for (i in 1:length(raw.txt)) raw.txt[i] = unlist(strsplit(raw.txt[i],"\t"))[2]
        
        # optional stage ("exclusivity pruning")
        #####################################################################################
        ####  get rid of abstracts that contain a keyword from other phenotypes/tissues/keywords
        if (exclusivity_pruning==TRUE) {
        
        # for pregnancy-related genes
        if (phe=="ENDOMETRIUM") regexp_not="myometr|([[:punct:]]|\\s)+cervi|([[:punct:]]|\\s)+uter[uaoi]+|placent"
        if (phe=="MYOMETRIUM") regexp_not="endometr|([[:punct:]]|\\s)+cervi|([[:punct:]]|\\s)+uter[uaoi]+|placent"
        if (phe=="UTERUS") regexp_not="endometr|myometr|([[:punct:]]|\\s)+cervi|placent"
        if (phe=="CERVIX") regexp_not="endometr|myometr|([[:punct:]]|\\s)+uter[uaoi]+|placent"
        if (phe=="PLACENTA") regexp_not="endometr|myometr|([[:punct:]]|\\s)+cervi|([[:punct:]]|\\s)+uter[uaoi]+"
        
        # for a control set of genes
        if (phe=="BLADDER") regexp_not="([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|([[:punct:]]|\\s)+dent|penile|prostat|trachea[ao]"
        if (phe=="BONE") regexp_not="bladder|([[:punct:]]|\\s)+dent|penile|prostat|trachea[ao]"
        if (phe=="DENTAL") regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|penile|prostat|trachea[ao]"
        if (phe=="PENILE") regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|([[:punct:]]|\\s)+dent|prostat|trachea[ao]"
        if (phe=="PROSTATE") regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|([[:punct:]]|\\s)+dent|penile|trachea[ao]"
        if (phe=="TRACHEA") regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|([[:punct:]]|\\s)+dent|penile|prostat"
        print(regexp_not)
        
        ## decide now whether ou want to use exclusivity filter or not
        bad = grep(regexp_not,raw.txt,ignore.case = T)  # ***
        raw.txt = raw.txt[-bad]    #                           ****
        }
        
        print(paste("number of abstacts (after exclusivity pruning): ",
                    ifelse(exclusivity_pruning,length(raw.txt),"NOT DONE"),sep=""))
        n_abstracts_2 = length(raw.txt)
        
        #####################################################################################
        ####  get rid of abstracts that contain other restricted code words
        bad.lines2=unique(grep("dementia|alzheim|schizo|adhd|autism|asperg",raw.txt,ignore.case=T))
        bad.lines3=unique(grep("purpura|fulminans|diabet|cancer|obesity|leukemi|fibroid",raw.txt,ignore.case=T))
        bad.lines4=unique(grep("chlamydia|coronar.{2,30}diseas| stroke|migraine|influenza",raw.txt,ignore.case=T))        
        bad.lines5=unique(grep("genetic defects|asthma|down syndrome|down's syndrome|klinefelt",raw.txt,ignore.case=T))
        bad.lines6=unique(grep("preimplantation genetic diagnosis|azoospermia",raw.txt,ignore.case=T))
        # not yet included:   tumor (but should not be used since "tumor necrosis factor".....)
        # not yet included: "syndrome" (too broad term?)
        
        bad.lines=unique(c(bad.lines2,bad.lines3,bad.lines4,bad.lines5,bad.lines6)); length(bad.lines)
        raw.txt=raw.txt[-bad.lines]
        print(paste("number of abstacts (after pop disease pruning): ",length(raw.txt),sep=""))
        n_abstracts_3 = length(raw.txt)
        
        #####################################################################################
        ####  get rid of abstracts that are not realistically important and might contain biases (ANIMAL studies)
        
        # identify bad abstracts
        n_abstracts = length(raw.txt)
        n_animals  = length(animals)
        lst = list()
        animal_test = NULL
        animal_regexp = paste("([[:punct:]]|\\s)+",animals,"([[:punct:]]|\\s)+",sep="")
        head(animal_regexp)
        for (i in 1:n_abstracts) {   # takes about a minute (per thousand abstracts)
                #print(paste(i,"/",n_abstracts,sep=""))
                tst = NULL
                for (j in 1:n_animals) {
                        tst=c(tst, length(grep(animal_regexp[j],raw.txt[i],ignore.case = T))>0)
                }
                animal_test=c(animal_test, sum(tst)>0)
                lst[[i]]=animals[which(tst)]
                rm(tst)
        }
        
        #table(animal_test)
        #sort(unique(unlist(lst)))
        
        # clean-up
        #raw.txt[which(animal_test)]
        raw.txt=raw.txt[-which(animal_test)]
        print(paste("number of abstacts (after animal pruning): ",length(raw.txt),sep=""))
        n_abstracts_4 = length(raw.txt)
        
        #####################################################################################
        ####  get rid of abstracts that contain forbiden medical terms (medical conditions)        
        
        n_abstracts = length(raw.txt)
        n_diseases  = length(diseases)
        lst = list()
        disease_test = NULL
        diseases_regexp = paste("([[:punct:]]|\\s)+",diseases,"([[:punct:]]|\\s)+",sep="")
        head(diseases_regexp)
        for (i in 1:n_abstracts) {   # takes about a minute (per thousand abstracts)
                #print(paste(i,"/",n_abstracts,sep=""))
                tst = NULL
                for (j in 1:n_diseases) {
                        
                        tst=c(tst, length(grep(diseases_regexp[j],raw.txt[i],ignore.case = T))>0)
                }
                disease_test=c(disease_test, sum(tst)>0)
                lst[[i]]=diseases[which(tst)]
        }
        
        #table(disease_test)
        #sort(unique(unlist(lst)))
        
        # final clean-up
        raw.txt=raw.txt[-which(disease_test)]
        n_abstracts_5 = length(raw.txt)
        print(paste("number of abstacts (after rare disease pruning): ",length(raw.txt),sep=""))
        # note that  "retractions" are taken care of in previous text mining script (by Julius)
        
        # DECISION WHETHER TRANSLATOR should be used
        if (translator_usage==TRUE) {
                translator = read.table("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/TRANSLATOR_misspelled_gene_names.txt",
                                        stringsAsFactors=F,h=T,sep="\t") # TRANSLATOR OF SOME GENE NAMES
                from_length = nchar(translator$from)
                translator = translator[order(from_length,decreasing = T),]
        } else { translator = data.frame(from="111111",to="222222") }

        
        # define dangerous gene names that are also biomed acronyms
        restricted_acronyms = unique(c( "AGA","SGA","LGA","FGR","AFD","ART","DMP","GO","IUGR","FGR","MC","DC","ECM","CDG",
                                        "SPTB","PTL","WAS","FTL","BPD","RDS","PTD","PTB","PROM","PPROM","PTL","CPHD","IAI",
                                        "NDN","TSL","POR","CAT","RAT","PIG","CATS", "MSC","PAH","PLEC","PIH","IVF","HRT",
                                        "CI","OR","RR","CC","SNP","MDR","RNA","DNA","ISCI","LOD","CAD","PGP","ROC","CPE",
                                        "MRI","CSM","HIV","HPV", "SDS","PAGE","SAGE", "FIGO", "ADO","PCR","QPCR","IVH","OI",
                                        "ROP","OS","RDS","BPD","ROP","AIM","THE", "PGD","ADO","SDS","PLEC","HUVEC","ERA",
                                        "SPARC","FOR","THE","BCM","HEEC","MSC","LNG","AMP","CERTL","DDT","ANOVA","COCP","MIAC",
                                        "BAD","PRL","PGF","TERT","CAC","CTC","TTC","ISH","ECS","ESC","MPA","CGB","CGA","EVT"))
        # almost all above mentioned acronyms are included in gene-name TRANSLATOR file and thus can be detected via their "long-name"
        #hg[which(hg$HUGO %in% restricted_acronyms),]
        # congenital disorder of glycosylation (CDG) #IAI = intraamniotic infection # OI = Osteogenesis Imperfecta
        
cumm1=NULL  # cummulation of potential gene names extracted from abstracts without TRANSLATOR
cumm2=NULL  # cummulation of REAL gene names extracted from abstracts using TRANSLATOR
cumm3=NULL  # cummulation of ALL unique-in-one-abstract gene names (ddtected with or without translator). for "times-mentioned" threshold


### perform cycling through abstracts with gene-name extraction procedure
for (j in 1:length(raw.txt)) { # takes a while...
        print(paste(j,"/",length(raw.txt),sep=""))
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
        
        # extract all remaining POTENTIAL gene names
        lst1=unlist(regmatches(txt, gregexpr('([[:punct:]]|[[:upper:]]|[0-9])+', txt)))
        lst2=unlist(regmatches(lst1, gregexpr('^[[:upper:]]{2,20}.*', lst1)))
        lst2 = sort(unique(lst2))
        lst2 = lst2[which(! lst2 %in% "REPLACED")]  # do not include artifacts

        # ELIMINATE GENES THAT ARE EASILY CONFUSED WITH THE TERMS FROM OBSTETRIC/MEDICAT/TECHNICAL FIELDS
        
        # eliminate from the first set
        lst3 = lst2[which(! lst2 %in% restricted_acronyms)]
        
        cumm1=c(cumm1,unique(lst3))
        cumm2=c(cumm2,unique(replaced))
        cumm3=c(cumm3,unique(c(lst3,replaced))) # will be used to estimate number of instances when a gene was mentioned
        rm(txt,lst1,lst2,lst3,replaced)
}


# directly extracted POTENTIAL genes
genes_1 = cumm1[which(cumm1 %in% hg_genes)]
hist(as.numeric(table(genes_1)),breaks=100,col="grey")
#genes_1 = sort(unique(genes_1))

# detected via TRANSLATOR genes
genes_2 = cumm2[which(cumm2 %in% hg_genes)]
hist(as.numeric(table(genes_2)),breaks=100,col="grey")
#genes_2 = sort(unique(genes_2))

# detected via BOTH METHODS
genes_3 = cumm3[which(cumm3 %in% hg_genes)]
hist(as.numeric(table(genes_3)),breaks=100,col="grey")
#genes_3 = sort(unique(genes_3))

library(gplots)
temp = list(extracted = genes_1,translated = genes_2)
venn(temp); rm(temp)

# all together
#if (length(genes_2)==0) {genes = sort(genes_1)} else {genes=sort(c(genes_1,genes_2))}
genes = genes_3

temp1=data.frame(gene=names(table(genes)),freq=as.numeric(table(genes)),stringsAsFactors = F)
temp2= merge(temp1,hg,by.x="gene",by.y="HUGO",all.x=T); rm(temp1)
gene.freq=temp2[rev(order(temp2$freq)),]; rm(temp2)
#gene.freq[1:20,]
gene_lists[[phe]]=gene.freq
rm(gene.freq, genes, genes_1,genes_2, cumm1,cumm2,cumm3)

summary_lists[[phe]] = c(n_abstracts_0,n_abstracts_1,n_abstracts_2,n_abstracts_3,n_abstracts_4,n_abstracts_5)
rm(n_abstracts_0,n_abstracts_1,n_abstracts_2,n_abstracts_3,n_abstracts_4,n_abstracts_5)

} # end of cycling through various phenotypes





############################################
################   save the results

# for pregnancy related genes

if( (exclusivity_pruning==TRUE)&(translator_usage==TRUE)) {
obg_xcl_trn = gene_lists  # obg = OBGYN, xcl = exclusivity filter ON, trn = TRANSLATOR ON
obg_xcl_trn_stats = summary_lists  # obg = OBGYN, xcl = exclusivity filter ON, trn = TRANSLATOR ON
obg_xcl_trn_hash = system(paste("git log --pretty=format:'%h' -n 1"),intern=TRUE)
}

if( (exclusivity_pruning==FALSE)&(translator_usage==TRUE)) {
obg_nxc_trn = gene_lists  # obg = OBGYN, nxc = exclusivity filter OFF, trn = TRANSLATOR ON
obg_nxc_trn_stats = summary_lists  # obg = OBGYN, nxc = exclusivity filter OFF, trn = TRANSLATOR ON
obg_nxc_trn_hash = system(paste("git log --pretty=format:'%h' -n 1"),intern=TRUE)
}

if( (exclusivity_pruning==TRUE)&(translator_usage==FALSE)) {
obg_xcl_unt = gene_lists  # obg = OBGYN, xcl = exclusivity filter ON, unt = TRANSLATOR OFF
obg_xcl_unt_stats = summary_lists  # obg = OBGYN, xcl = exclusivity filter ON, unt = TRANSLATOR OFF
obg_xcl_unt_hash = system(paste("git log --pretty=format:'%h' -n 1"),intern=TRUE)
}

if( (exclusivity_pruning==FALSE)&(translator_usage==FALSE)) {
obg_nxc_unt = gene_lists  # obg = OBGYN, nxc = exclusivity filter OFF, unt = TRANSLATOR OFF
obg_nxc_unt_stats = summary_lists  # obg = OBGYN, nxc = exclusivity filter OFF, unt = TRANSLATOR OFF
obg_nxc_unt_hash = system(paste("git log --pretty=format:'%h' -n 1"),intern=TRUE)
}

# short version
save(list=c("obg_xcl_trn","obg_nxc_trn",
            "obg_xcl_trn_stats","obg_nxc_trn_stats",
            "obg_xcl_trn_hash","obg_nxc_trn_hash"),
     file="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/obgn_genes.RData")

# long version
save(list=c("obg_xcl_trn","obg_xcl_unt","obg_nxc_trn","obg_nxc_unt",
            "obg_xcl_trn_stats","obg_xcl_unt_stats","obg_nxc_trn_stats","obg_nxc_unt_stats",
            "obg_xcl_trn_hash","obg_xcl_unt_hash","obg_nxc_trn_hash","obg_nxc_unt_hash"),
     file="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/obgn_genes.RData")


# for control set of genes
if( (exclusivity_pruning==TRUE)&(translator_usage==TRUE)) {
        ctrl_xcl_trn = gene_lists  # ctrl = CONTROL, xcl = exclusivity filter ON, trn = TRANSLATOR ON
        ctrl_xcl_trn_stats = summary_lists  # ctrl = CONTROL, xcl = exclusivity filter ON, trn = TRANSLATOR ON
        ctrl_xcl_trn_hash = system(paste("git log --pretty=format:'%h' -n 1"),intern=TRUE)
}

if( (exclusivity_pruning==FALSE)&(translator_usage==TRUE)) {
        ctrl_nxc_trn = gene_lists  # ctrl = CONTROL, nxc = exclusivity filter OFF, trn = TRANSLATOR ON
        ctrl_nxc_trn_stats = summary_lists  # ctrl = CONTROL, nxc = exclusivity filter OFF, trn = TRANSLATOR ON
        ctrl_nxc_trn_hash = system(paste("git log --pretty=format:'%h' -n 1"),intern=TRUE)
}


save(list=c("ctrl_xcl_trn","ctrl_nxc_trn",
            "ctrl_xcl_trn_stats","ctrl_nxc_trn_stats",
            "ctrl_xcl_trn_hash","ctrl_nxc_trn_hash"),
     file="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/cntrl_genes.RData")




############################################
###########   PREVIEW

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


########################################################################
########################################################################
######################### export the gene sets

rm(list=ls())  # cleanup
load("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/obgn_genes.RData")
#load("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/WORK_FILES/cntrl_genes.RData")
lst=ls()
obj_lst = lst[grep("^obg_.{7,7}$|^ctrl_.{7,7}$",lst)]  # gene frequencies per each phenotype/tissue

for (z in 1:length(obj_lst)) {  # for each type of settings
        temp_obj = get(obj_lst[z])
        types = names(temp_obj)
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






################################################################
################################################################
################################################################
################################################################

......      not finished below   .........

#############   EXTARCT OVERLAPING AND UNIQUE GENES for each phenotype
all_genes = list()
individual = list()
phenos =  c("CERVIX","ENDOMETRIUM","MYOMETRIUM","UTERINE","PLACENTA")
for (pheno in phenos) {
#pheno = "CERVIX"
phes = phenos[phenos != pheno]
lll = NULL; for (phe in phes) lll=c(lll, unlist(gene_lists[[phe]]$gene))
ooo = gene_lists[[pheno]]$gene
individual[[pheno]] = ooo[which( ! ooo %in% sort(unique(lll)))]
all_genes[[pheno]] = unique(unlist(ooo))
rm(lll,ooo,phes)
}
all_unq =  unique(unlist(individual)) 
all_genes = unique(unlist(all_genes))
common = unique(all_genes[ ! all_genes %in% all_unq ])


new_collection=NULL
for (i in 1:(length(phenos)+1))  {
        if (i<=length(phenos)) {
        #        i=1
        gene_names = data.frame(V1 = unique(individual[[phenos[i]]]))
        m = merge(hg , gene_names, by.x= "HUGO" ,by.y= "V1", all= F)
        temp = data.frame(ENTREZ= m$ENTREZ, GENESET = paste("PubMed:",phenos[i],sep=""),Descript=".")
        new_collection = rbind(new_collection,temp)
        rm(temp,m,gene_names)
        } else {
                gene_names = data.frame(V1 = unique(common))
                m = merge(hg , gene_names, by.x= "HUGO" ,by.y= "V1", all= F)
                temp = data.frame(ENTREZ= m$ENTREZ, GENESET = "PubMed:common",Descript=".")
                new_collection = rbind(new_collection,temp)
                rm(temp,m,gene_names)
        }
}

dim(new_collection); head(new_collection)
table(new_collection$GENESET)
table(collection$GENESET)
new_printout=new_collection
colnames(new_printout)[1]=paste("##",colnames(new_printout)[1],sep="")
write.table(new_printout,"~/Biostuff/MOBA_GESTAGE_GWAS/INRICH/PubMed_OBGYN_5tissuesUnqComm.txt",
            row.names=F,col.names=T,sep="\t",quote=F)

table(new_collection$GENESET)


my_genes = c("SP3", "SLC1A1", "TLR4", "IGF2", "MMP9" ,"EDN3")
colnames(new_collection)[1]="ENTREZ"
ggg=data.frame(V1=new_collection$ENTREZ[which(new_collection$GENESET=="PubMed:ENDOMETRIUM")])
mmm = merge(hg , ggg, by.x= "ENTREZ" ,by.y= "V1", all= F)
my_genes %in% mmm$HUGO
my_genes %in% common



