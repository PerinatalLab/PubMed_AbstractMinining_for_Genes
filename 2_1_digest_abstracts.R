#!/usr/bin/Rscript

# description: extracts/mines PubMed abstracts that are related to PREGNANCY 
# but are not biased by animal studies and are +/- unique to the specific tissue


# 1) search the PubMed website with code words **tissue**  AND **gene**
# 2) download the text file
# 3) run Julius' script that that awk-eliminates everything except abstracts
# 4) run this script that further prunes/cleanes abstracts 


###  get typical disease-names from definitions in ICD-code-file 
#    .. (in order to get rid of abstracts with other diseases)
txt = read.table("CMS32_DESC_SHORT_DX.txt",h=F,sep="\t",stringsAsFactors = F) #~/Biostuff/ICD-9-CM-v32-master-descriptions/
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

##  get the list of animal-related terms that should not be in the abstracts (non-human subjects)
anim1 = read.table("animal_name_indicators.txt",stringsAsFactors = F,h=F,sep="\t") #~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY/
anim1=sort(unique(tolower(anim1[,1])))
anim2 = unique(c("animal","cattle","buffalo","ruminant", "cow","dog","rat","pig","cat","lion","horse","monkey","mouse",
                 "cows","dogs","rats","pigs","cats","lions","horses","monkeys","mice"))
animals=sort(unique(c(anim1,anim2))); rm(anim1,anim2)


PubMedDir="./PubMed_DIGEST/" #~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY
pttrn = "PLACEN|CERVIX|MYOMETR|ENDOMETR|UTER|PENILE|BLADD|BONE|PROSTAT|SKIN|INTESTIN|MUSCLE|LIVER|LUNG"
file_list = list.files(PubMedDir,pattern=pttrn)
files_ok = file_list[grep("abstracts",file_list)]
pheno=NULL; for (i in 1:length(files_ok))pheno = c(pheno, unlist(strsplit(files_ok[i],"_"))[4]); print(pheno)

phes = pheno
cleaned_abstracts_exclusivityON = list() # cummulator of cleaned abstracts
cleaned_abstracts_exclusivityOFF = list() # cummulator of cleaned abstracts
stats_excl_ON = stats_excl_OFF = NULL # cummulator of cleaning summary stats

for (exclusivity_pruning in c(TRUE,FALSE)) {
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
                        if (phe=="ENDOMETRIUM") regexp_not="myometr|([[:punct:]]|\\s)+cervi|([[:punct:]]|\\s)+uter[uaoi]+"
                        if (phe=="MYOMETRIUM") regexp_not="endometr|([[:punct:]]|\\s)+cervi|([[:punct:]]|\\s)+uter[uaoi]+"
                        if (phe=="UTERUS") regexp_not="endometr|myometr|([[:punct:]]|\\s)+cervi"
                        if (phe=="CERVIX") regexp_not="endometr|myometr|([[:punct:]]|\\s)+uter[uaoi]+"
                        if (phe=="PLACENTA") regexp_not="endometr|myometr|([[:punct:]]|\\s)+cervi|([[:punct:]]|\\s)+uter[uaoi]+"
                        
                        # for a control set of genes
                        if (phe=="BLADDER")     regexp_not="([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|penile|prostat|intestin|([[:punct:]]|\\s)+liver|([[:punct:]]|\\s)+hepat|([[:punct:]]|\\s)+lung|pulmonary|muscle|([[:punct:]]|\\s)+skin |([[:punct:]]|\\s)+derma"
                        if (phe=="BONE")        regexp_not="bladder|penile|prostat|intestin|([[:punct:]]|\\s)+liver|([[:punct:]]|\\s)+hepat|([[:punct:]]|\\s)+lung|pulmonary|muscle|([[:punct:]]|\\s)+skin |([[:punct:]]|\\s)+derma"
                        if (phe=="PENILE")      regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|prostat|intestin|([[:punct:]]|\\s)+liver|([[:punct:]]|\\s)+hepat|([[:punct:]]|\\s)+lung|pulmonary|muscle|([[:punct:]]|\\s)+skin |([[:punct:]]|\\s)+derma"
                        if (phe=="PROSTATE")    regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|penile|intestin|([[:punct:]]|\\s)+liver|([[:punct:]]|\\s)+hepat|([[:punct:]]|\\s)+lung|pulmonary|muscle|([[:punct:]]|\\s)+skin |([[:punct:]]|\\s)+derma"
                        if (phe=="INTESTINE")   regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|penile|prostat|([[:punct:]]|\\s)+liver|([[:punct:]]|\\s)+hepat|([[:punct:]]|\\s)+lung|pulmonary|muscle|([[:punct:]]|\\s)+skin |([[:punct:]]|\\s)+derma"
                        if (phe=="LIVER")       regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|penile|prostat|intestin|([[:punct:]]|\\s)+lung|pulmonary|muscle|([[:punct:]]|\\s)+skin |([[:punct:]]|\\s)+derma"
                        if (phe=="LUNGS")       regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|penile|prostat|intestin|([[:punct:]]|\\s)+liver|([[:punct:]]|\\s)+hepat|muscle|([[:punct:]]|\\s)+skin |([[:punct:]]|\\s)+derma"
                        if (phe=="MUSCLE")      regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|penile|prostat|intestin|([[:punct:]]|\\s)+liver|([[:punct:]]|\\s)+hepat|([[:punct:]]|\\s)+lung|pulmonary|([[:punct:]]|\\s)+skin |([[:punct:]]|\\s)+derma"
                        if (phe=="SKIN")        regexp_not="bladder|([[:punct:]]|\\s)+oste[oa]|([[:punct:]]|\\s)+bone|penile|prostat|intestin|([[:punct:]]|\\s)+liver|([[:punct:]]|\\s)+hepat|([[:punct:]]|\\s)+lung|pulmonary|muscle"
                        
                        print(regexp_not)
                        
                        ## decide now whether ou want to use exclusivity filter or not
                        bad = grep(regexp_not,raw.txt,ignore.case = T)
                        raw.txt = raw.txt[-bad]    
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
                
                # clean-up
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
              
                # final clean-up
                raw.txt=raw.txt[-which(disease_test)]
                n_abstracts_5 = length(raw.txt)
                print(paste("number of abstacts (after rare disease pruning): ",length(raw.txt),sep=""))
                # note that  "retractions" are taken care of in previous text mining script (by Julius)
                
                # saving
                if(exclusivity_pruning == TRUE) {
                        cleaned_abstracts_exclusivityON[[phe]] = raw.txt
                        nmbrs = data.frame(original=n_abstracts_0, longtxt=n_abstracts_1,
                                           exclusivity=n_abstracts_2, popDisease=n_abstracts_3,
                                           animals=n_abstracts_4, rareDisease=n_abstracts_5,row.names=phe)
                        stats_excl_ON = rbind(stats_excl_ON,nmbrs); rm(nmbrs)
                }
                
                # saving
                if(exclusivity_pruning == FALSE) {
                        cleaned_abstracts_exclusivityOFF[[phe]] = raw.txt
                        nmbrs = data.frame(original=n_abstracts_0, longtxt=n_abstracts_1,
                                           exclusivity=n_abstracts_2, popDisease=n_abstracts_3,
                                           animals=n_abstracts_4, rareDisease=n_abstracts_5,row.names=phe)
                        stats_excl_OFF = rbind(stats_excl_OFF,nmbrs); rm(nmbrs)
                }
                
                rm(raw.txt)
                
        }
}


cleaned_abstracts_exclusivityON[["stats"]] = stats_excl_ON
cleaned_abstracts_exclusivityOFF[["stats"]] = stats_excl_OFF

#  save what was generated (cleaned)
out_dir="./PubMed_PRUNE/" #~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY
save(list=c("cleaned_abstracts_exclusivityON","cleaned_abstracts_exclusivityOFF"),
     file=paste(out_dir,"cleaned_abstracts_GYNECOLOGY.RData",sep=""))

rm(list=ls())

# load what was generated (cleaned)
#out_dir="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun_GYNECOLOGY/PubMed_PRUNE/"
#load(paste(out_dir,"cleaned_abstracts_GYNECOLOGY.RData",sep=""))

