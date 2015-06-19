
# extract PubMed genes that are related to PREGNANCY (based on abstract mining)

# 1) search the PubMed website with code words **tissue**  AND **gene**, download the text file
# 2) run Julius' script that that eliminates everything except abstracts
# 3) run this script that further prunes abstracts and extracts gene names


# load the full collection human genes
hg=read.table("~/Biostuff/hg19_HUMAN_GENES/ucsc_HUGO_ENTREZ_chr1-23_withDescriptions_hg19_PROCESSED.txt",
              stringsAsFactors=F,h=F); dim(hg); head(hg); table(hg$V1)
colnames(hg)=c("CHR","START","END","ENTREZ","HUGO","Description")
hg_genes=hg[which(nchar(hg$HUGO)>2),"HUGO"]; length(hg_genes) # arguably too much risk with short gene names
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
animals=sort(unique(c(anim1,anim2)))


PubMedDir="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_DIGEST/"
file_list = list.files(PubMedDir,pattern="PLACEN|CERVIX|MYOMETR|ENDOMETR|UTER")  # pregnancy-related genes
#file_list = list.files(PubMedDir,pattern="PENILE|BLADD|BONE|DENTAL|PROSTAT|TRACHE")  # control set of genes (other tissues)
files_ok = file_list[grep("abstracts",file_list)]
pheno=NULL; for (i in 1:length(files_ok))pheno = c(pheno, unlist(strsplit(files_ok[i],"_"))[4]); print(pheno)

#phes = c("PREGNANCY","CERVIX","ENDOMETRIUM","MYOMETRIUM","PRETERM","BIRTH","GESTATIONAL","UTERINE","PLACENTA")
#phes =  c("ENDOMETRIUM","MYOMETRIUM","CERVIX","UTERUS","PLACENTA")
phes = pheno


##  number of Abstracts before the cleaning
n_abs0 = NULL
for (phe in phes) {
        #phe = "ENDOMETRIUM"
        file_name = files_ok[which(pheno==phe)]
        raw.txt=readLines(paste(PubMedDir,file_name,sep="")) ; length(raw.txt)
        #raw.txt= raw.txt[ grep("^Abstract",raw.txt)]  # only use content of the Abstract (due to exclusions of term words)
        n_abs0 = c(n_abs0,length(raw.txt))
}

        
gene_lists = list()
collection1= NULL  # for those with no limit on gene frequency (number of abstracts containing that gene name)
collection2= NULL  # for those with minimum 2 abstracts mentioning the gene name

for (phe in phes) {
        print(phe)
        file_name = files_ok[which(pheno==phe)]
        raw.txt=readLines(paste(PubMedDir,file_name,sep=""))
        print(paste("number of abstacts (initial): ",length(raw.txt),sep=""))
        
        goo_length = which(nchar(raw.txt)>100)
        raw.txt=raw.txt[goo_length]; rm(goo_length)
        print(paste("number of abstacts (length > 100 smbls): ",length(raw.txt),sep=""))
        
        # get rid of tab symbol in the begining of the text string
        for (i in 1:length(raw.txt)) raw.txt[i] = unlist(strsplit(raw.txt[i],"\t"))[2]
        
        # optional stage:
        #####################################################################################
        ####  get rid of abstracts that contain a keyword from other phenotypes/tissues/keywords
        
        if (phe=="ENDOMETRIUM") regexp_not="myometr|([[:punct:]]|\\s)+cervi|([[:punct:]]|\\s)+uter[uaoi]+|placent"
        if (phe=="MYOMETRIUM") regexp_not="endometr|([[:punct:]]|\\s)+cervi|([[:punct:]]|\\s)+uter[uaoi]+|placent"
        if (phe=="UTERUS") regexp_not="endometr|myometr|([[:punct:]]|\\s)+cervi|placent"
        if (phe=="CERVIX") regexp_not="endometr|myometr|([[:punct:]]|\\s)+uter[uaoi]+|placent"
        if (phe=="PLACENTA") regexp_not="endometr|myometr|([[:punct:]]|\\s)+cervi|([[:punct:]]|\\s)+uter[uaoi]+"
        print(regexp_not)
        
        bad = grep(regexp_not,raw.txt)
        raw.txt = raw.txt[-bad]
        print(paste("number of abstacts (after exclusivity pruning): ",length(raw.txt),sep=""))
        
        #####################################################################################
        ####  get rid of abstracts that contain other restricted code words
        
        bad.lines3=unique(grep("purpura|fulminans|diabet|cancer|obesity|leukemi|alzheim|schizo|adhd|fibroid",raw.txt,ignore.case=T))
        bad.lines4=unique(grep("chlamydia|coronar.{2,30}diseas| stroke|migraine",raw.txt,ignore.case=T))        
        bad.lines5=unique(grep("genetic defects|asthma|down syndrome",raw.txt,ignore.case=T))
        bad.lines6=unique(grep("preimplantation genetic diagnosis|azoospermia",raw.txt,ignore.case=T))
        # not yet included:   tumor (but should not be used since "tumor necrosis factor".....)
        
        bad.lines=unique(c(bad.lines3,bad.lines4,bad.lines5,bad.lines6)); length(bad.lines)
        raw.txt=raw.txt[-bad.lines]
        print(paste("number of abstacts (after pop disease pruning): ",length(raw.txt),sep=""))
        
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
        
        # clean-up
        raw.txt=raw.txt[-which(disease_test)]
        print(paste("number of abstacts (after rare disease pruning): ",length(raw.txt),sep=""))
        
        # note that  "retractions" are taken care of in previous text mining script (by Julius)

        # TRANSLATOR OF SOME GENE NAMES
        translator= read.table("~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/TRANSLATOR_misspelled_gene_names.txt",
                               stringsAsFactors=F,h=T,sep="\t")
        #head(translator)
        from_length = nchar(translator$from)
        translator = translator[order(from_length,decreasing = T),]
        #head(translator); rm(from_length)

        # IF YOU DO NOT WANT TO USE TRANSLATOR - activate the following line:
        #translator = data.frame(from="111111",to="222222")

cumm1=NULL  # cummulation of potential gene names extracted from abstracts
cumm2=NULL  # cummulation of REAL gene names extracted from abstracts using TRANSLATOR

for (j in 1:length(raw.txt)) {
        
        # this is the text that we will be working with in this cycle
        txt=raw.txt[j]
        
        ##  use the TRANSLATOR to EXTARCT (!) gene names that are written in a non-standard manner
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

        # remove all possible punctuation marks, bt preserve the hyphen
        txt=paste(unlist(strsplit(txt,"-")),collapse="zzzzz") # "zubiquitilation of hyphens"
        txt=paste(unlist(strsplit(txt,"[[:punct:]]")),collapse=" ")
        txt=paste(unlist(strsplit(txt,"zzzzz")),collapse="-")
        
        # extract all remaining POTENTIAL gene names
        lst1=unlist(regmatches(txt, gregexpr('([[:punct:]]|[[:upper:]]|[0-9])+', txt)))
        lst1=unlist(regmatches(lst1, gregexpr('^[[:upper:]]{2,20}.*', lst1)))
        lst1 = sort(unique(lst1))
        lst1 = lst1[which(!lst1 %in% "REPLACED")]

        cumm1=c(cumm1,unique(lst1))
        cumm2=c(cumm2,unique(replaced))
        rm(txt,lst1,replaced)
}


# directly extracted POTENTIAL genes
genes_1 = cumm1[which(cumm1 %in% hg_genes)]
hist(as.numeric(table(genes_1)),breaks=100,col="grey")
#genes_1 = sort(unique(genes_1))

# detected via TRANSLATOR genes
genes_2 = cumm2[which(cumm2 %in% hg_genes)]
#hist(as.numeric(table(genes_2)),breaks=100,col="grey")
#genes_2 = sort(unique(genes_2))

library(gplots)
temp = list(extracted = genes_1,translated = genes_2)
venn(temp); rm(temp)

# ELIMINATE GENES (from the genes_1 set) THAT ARE EASILY CONFUSED WITH THE TERMS FROM OBSTETRIC FIELD
# define dangerous gene names that are also biomed acronyms
restricted_acronyms = unique(c( "AGA","SGA","LGA","FGR","AFD","ART","DMP","GO","IUGR","FGR","MC","DC","ECM",
                                "SPTB","PTL","WAS","FTL","BPD","RDS","PTD","PTB","PROM","PPROM","PTL","CPHD",
                                "NDN","TSL","POR","CAT","RAT","PIG","CATS", "MSC","PAH","PLEC","PIH","IVF","HRT",
                                "CI","OR","RR","CC","SNP","MDR","RNA","DNA","ISCI","LOD","CAD","PGP","ROC","CPE",
                                "MRI","CSM","HIV","HPV", "SDS","PAGE","SAGE", "FIGO", "ADO","PCR","QPCR","IVH",
                                "ROP","OS","RDS","BPD","ROP","AIM","THE", "PGD","ADO","SDS","PLEC","HUVEC","ERA",
                                "SPARC","FOR","THE","BCM","HEEC","MSC","LNG","AMP","CERTL","DDT","ANOVA","COCP",
                                "BAD","PRL","PGF","TERT","CAC","CTC","TTC","ISH","ECS","ESC","MPA","CGB","CGA","EVT"))
# congenital disorder of glycosylation (CDG)
#IAI - intraamniotic infection
#Osteogenesis Imperfecta
#MIAC

# eliminate from the first set
genes_1 = genes_1[which( ! genes_1 %in% restricted_acronyms)]

# all together
if (length(genes_2)==0) {genes = sort(genes_1)} else {genes=sort(c(genes_1,genes_2))}

temp1=data.frame(gene=names(table(genes)),freq=as.numeric(table(genes)),stringsAsFactors = F)
temp2= merge(temp1,hg,by.x="gene",by.y="HUGO",all.x=T); rm(temp1)
gene.freq=temp2[rev(order(temp2$freq)),]; rm(temp2)
#gene.freq[1:20,]

#dim(gene.freq)
#hist(as.numeric(table(genes)),breaks=100,col="grey")

p1= data.frame(ENTREZ = gene.freq[ gene.freq$freq>=1,"ENTREZ"], GENESET = paste("PM",":",phe,sep=""),Descript="obs>=1")
p2= data.frame(ENTREZ = gene.freq[ gene.freq$freq>=2,"ENTREZ"], GENESET = paste("PM",":",phe,sep=""),Descript="obs>=2")
#p3= data.frame(ENTREZ = gene.freq[ gene.freq$freq>=3,"ENTREZ"], GENESET = paste("PubMed",":",phe,"_3plus",sep=""),Descript="Abstract mining")
#p4= data.frame(ENTREZ = gene.freq[ gene.freq$freq>=4,"ENTREZ"], GENESET = paste("PubMed",":",phe,"_4plus",sep=""),Descript="Abstract mining")
#temp = rbind(p1,p2,p3,p4)
#head(p1); head(collection)
collection1 = rbind( collection1, p1) 
collection2 = rbind( collection2, p2) 
#dim(collection); head(collection)


#####

#list.files(".",pattern = "^HAR")
gene_lists[[phe]]=gene.freq
rm(gene.freq, genes, genes_1,genes_2, cumm1,cumm2,p1,p2)
} # end of cycling through various phenotypes




###########   PREVIEW
temp=list(endm=gene_lists$ENDOMETRIUM$gene,miom=gene_lists$MYOMETRIUM$gene,
          cerv = gene_lists$CERVIX$gene,uter = gene_lists$UTERUS$gene,
          plac = gene_lists$PLACENTA$gene)
venn(temp)


# for control-set of tissues/phenotypes
temp=list(bladd=gene_lists$BLADDER$gene,bone=gene_lists$BONE$gene,
          penile = gene_lists$PENILE$gene,prost = gene_lists$PROSTATE$gene,
          trach = gene_lists$TRACHEA$gene)
venn(temp)



da = gene_lists$TRACHEA
head(da)


###############   create a simplistic gene-set file  ( where sets might overlap)
table(collection1$GENESET)
table(collection2$GENESET)

###   save all FIVE tissues (including fetal placenta)
# no restriction on gene frequency
printout=collection1
colnames(printout)[1]=paste("##",colnames(printout)[1],sep="")
folder = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_GENES/"
write.table(printout,paste(folder,"PubMed_OBGuntrn1min_5tissOverlap.txt",sep=""),
            row.names=F,col.names=T,sep="\t",quote=F) # translated or untranslated ? 

# with restriction on gene frequency
printout=collection2
colnames(printout)[1]=paste("##",colnames(printout)[1],sep="")
folder = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_GENES/"
write.table(printout,paste(folder,"PubMed_OBGuntrn2min_5tissOverlap.txt",sep=""),
            row.names=F,col.names=T,sep="\t",quote=F) # translated or untranslated ? 


###   now without PLACENTA ( since that is a fetal tissue)
# no restriction on gene frequency
printout=collection1[which(collection1$GENESET != "PM:PLACENTA"),]
colnames(printout)[1]=paste("##",colnames(printout)[1],sep="")
folder = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_GENES/"
write.table(printout,paste(folder,"PubMed_OBGuntrn1min_4tissOverlap.txt",sep=""),
            row.names=F,col.names=T,sep="\t",quote=F) # translated or untranslated ? 

# with a restriction on gene frequency
printout=collection2[which(collection2$GENESET != "PM:PLACENTA"),]
colnames(printout)[1]=paste("##",colnames(printout)[1],sep="")
folder = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_GENES/"
write.table(printout,paste(folder,"PubMed_OBGuntrn2min_4tissOverlap.txt",sep=""),
            row.names=F,col.names=T,sep="\t",quote=F) # translated or untranslated ? 


### save the control-set of genes
# no restriction on gene frequency
printout=collection1
colnames(printout)[1]=paste("##",colnames(printout)[1],sep="")
folder = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_GENES/"
write.table(printout,paste(folder,"PubMed_CONTROLtrnsl1min_6tissOverlap.txt",sep=""),
            row.names=F,col.names=T,sep="\t",quote=F) # translated or untranslated ? 

# with restriction on gene frequency
printout=collection2
colnames(printout)[1]=paste("##",colnames(printout)[1],sep="")
folder = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_GENES/"
write.table(printout,paste(folder,"PubMed_CONTROLtrnsl2min_6tissOverlap.txt",sep=""),
            row.names=F,col.names=T,sep="\t",quote=F) # translated or untranslated ? 





################################################################
################################################################
################################################################
################################################################

head(collection)

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




##### raw texts

raw.txt[171:185]


