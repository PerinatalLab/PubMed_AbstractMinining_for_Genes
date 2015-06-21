

###     this script (semi)automatically expands the TRANSLATOR file
#       by detecting gene names that are often erroneously abbreviated 
#       (by Jonas 2015 June 20)


###  get as many abstracts as possible (for later use)
PubMedDir="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/PubMed_DIGEST/"
file_list = list.files(PubMedDir,pattern="abstracts")
raw.txt = NULL
for (n in 1:length(file_list)) {
        tmp=readLines(paste(PubMedDir,file_list[n],sep=""))
        raw.txt=c(raw.txt,tmp); rm(tmp)
}
raw.txt=raw.txt[which(nchar(raw.txt)>100)]; length(raw.txt)


# get the full list of human genes (for immediate use)
hg=read.table("~/Biostuff/hg19_HUMAN_GENES/ucsc_HUGO_ENTREZ_chr1-23_withDescriptions_hg19_PROCESSED.txt",
              stringsAsFactors=F,h=F); dim(hg); head(hg); table(hg$V1)
colnames(hg)=c("CHR","START","END","ENTREZ","HUGO","Description")


## get the frequency of known gene names that start with 2/3/4 letters, followed by digits (RUN it separately for 2/3/4 letters ***)
ixs=grep("^[A-Z]{4,4}[0-9]+",hg$HUGO) # grab 2/3/4-letter prefixes ***
txt=hg[ixs,"HUGO"]
txt=txt[which(nchar(txt)<=(4+3))] # we are only interesting in genes that are like XYZW123
lst2=unlist(regmatches(txt, gregexpr("^[A-Z]{4,4}", txt)))  # grab 2/3/4-letter prefixes ***
words = names(table(lst2)); counts = as.numeric(table(lst2))
dat = data.frame(word=words,freq=counts,stringsAsFactors = F)
dat=dat[order(dat$freq,decreasing = T),]
dim(dat); dat[1:10,]


### THOROUGHLY PERFORM THE SEARCH OF GENE ROOT (alphabetic part) in all Abstracts

#report="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/geneTrnsltrFile_precursor_2letters.txt" # ***
#report="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/geneTrnsltrFile_precursor_3letters.txt" # ***
report="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/geneTrnsltrFile_precursor_4letters.txt"  # ***
for (rix in 1:dim(dat)[1]) {
        print(paste(rix,"/",dim(dat)[1],sep=""))
        offset=25  # how many characters to the left and right should be shown
        term_word = dat[rix,1] #"MMP"
        regexp_term=paste("([[:punct:]]|\\s)+",term_word,"([[:punct:]]|\\s)+",sep="")
        ixs = grep(regexp_term,raw.txt)
        txtout=NULL
        for (ix in ixs) {
                obj=gregexpr(regexp_term, raw.txt[ix])[[1]]
                from = ifelse(obj[1]<offset,0,obj[1]-offset)
                to = from+attr(obj,"match.length")+offset+offset
                txtout = c(txtout,substr(raw.txt[ix],from,to))
        }
        
        if (length(txtout)>10) {  # arbitrary threshold to claim a "grammar-gene"
                write(paste(dat[rix,],collapse="   _______________________   "),file=report,append=T)
                write(txtout[1:5],file=report,append=T)
                write(" ",file=report,append=T)
        } 
rm(txtout)
}



### NOT NECESSARY :  detailed preview of the extracted names and situations in which these names occur...
report="~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/geneTrnsltrFile_precursor_3letters.txt"
ttt=readLines(report)
selected = substr(ttt[grep("[A-Z]{3,3}.*____",ttt)],1,3) # change numbers ***
for (ix in 1:length(selected)) {
        offset=25
        term_word = selected[ix]
        regexp_term=paste("([[:punct:]]|\\s)+",term_word,"([[:punct:]]|\\s)+",sep="")
        zs = grep(regexp_term,raw.txt)
        txtout=NULL
        for (z in zs) {
                obj=gregexpr(regexp_term, raw.txt[z])[[1]]
                from = ifelse(obj[1]<offset,0,obj[1]-offset)
                to = from+attr(obj,"match.length")+offset+offset
                txtout = c(txtout,substr(raw.txt[z],from,to))
        }
        print(txtout)      
}



# define a function that creates a geneName TRANSLATOR table, using likely combinations of abbreviation strategies
trnsl_fun = function(roo,nrs) {  # roo = "FGF"; nrs = c(1,2,3)
        from1 = paste(roo,nrs,sep=" "); from2 = paste(roo,nrs,sep="-")
        from3 = paste("(",roo,")-",nrs,sep=""); from4 = paste("(",roo,") ",nrs,sep="")
        from5 = paste("(",roo,")",nrs,sep=""); to1_5 = rep(paste(roo,nrs,sep=""), 5)
        from1_5 = c(from1,from2,from3,from4,from5)
        piece = data.frame(from = from1_5, to = to1_5)
        return(piece)
}


#    these sets are created semiautomatically (using the above script mining Abstracts), with some ...
#  ... manual checking (to make sure that it is a real gene name)
roots2 = c("IL","SH","NF","NT")
roots3 = c("MMP","CDK","COX","BMP","TCF","BCL","NAT","TLR","IRF","MAP","ATF","CPT","EGR","IGF","NOS","PAR","PTH",
          "SOD","ADM","SDF","DMP","HAT","HSD","LAP","MRP","TSG","ADM") # created semi automatically, with manual checking
roots4 = c("SSTR","HDAC","HOXA","PARP","ADAM","GATA","ICAM","FGFR","SUMO","TIMP","RUNX","BRCA","DNMT","VCAM")
roots = c(roots2,roots3,roots4)

# this chunk cycles to all suspicious gene names and creates a translator table (only for real gene names with simple name)
TRANLSATOR = NULL
for (roo in roots) { # for each gene class
        tmp1 = sort(unique(hg[grep(paste("^",roo,"[0-9]{1,3}$",sep=""),hg$HUGO),"HUGO"])) # extract simple gene names ( = ABC230)
        frs = rep(nchar(roo)+1,length(tmp1)); tos = nchar(tmp1) # define start and end of numeric part
        nrs = NULL; for (r in 1:length(tmp1)) nrs = c(nrs, substr(tmp1[r],frs[r],tos[r])) # extract numberic part of the gene name
        nrs = sort(unique(as.numeric(nrs))); rm(tmp1,frs,tos) # cleanup
        if (length(nrs)>0) TRANLSATOR = rbind(TRANLSATOR, trnsl_fun(roo,nrs))
}


dim(TRANLSATOR)
TRANLSATOR[grep("MMP",TRANLSATOR$from),]

# save the newly generated translator table with the previouss translator table
translator.file = "~/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/TRANSLATOR_misspelled_gene_names.txt"
old_one = read.table(translator.file,h=T,sep="\t",stringsAsFactors = F)
new_one = rbind(old_one,TRANLSATOR)
good_lines = which(! duplicated(new_one$from))
new_one = new_one[good_lines,]
ord = nchar(new_one$from)
new_one = new_one[order(ord,decreasing = T),]
head(new_one); tail(new_one); dim(new_one)
write.table(new_one,translator.file,row.names=F,col.names=T,sep="\t",quote=F)

