#!/usr/bin/Rscript
# extract PubMed genes that are related to PREGNANCY (based on abstract mining)

# 1) search the PubMed website with code words **tissue**  AND **gene**
# 2) download the text file
# 3) run Julius' script that that eliminates everything except abstracts
# 4) run "2_1..." script that further prunes abstracts
# 5) run this script which extracts valid gene names

## USAGE: ./2_2_extract_gene_names.R path_to_dir_with_PubMed_subdirs path_to_restricted_acronym_file

args = commandArgs(TRUE)
working_dir = args[1]
acro_file = args[2]

# load the full collection of human genes
# CHANGE PATH HERE AS NEEDED
hg = read.table("ucsc_HUGO_ENTREZ_chr1-23_withDescriptions_hg19_PROCESSED.txt",
              stringsAsFactors=F, h=F)
colnames(hg)=c("CHR","START","END","ENTREZ","HUGO","Description")
# min length of gene name is now 2. might want to increase, to reduce risk of incorporating other acronyms
hg_genes=hg[which(nchar(hg$HUGO)>1),"HUGO"]; length(hg_genes)

# load dangerous gene names that are also biomed acronyms (PCR, DNA...)
acronyms_tbl = read.table(acro_file,sep="\t",h=T,stringsAsFactors=F)
restricted_acronyms = sort(unique(acronyms_tbl$Acronym))

# some acronyms are present in the gene-name TRANSLATOR file..
# ..and thus can be detected via their "long-name"

# load what was generated in previous script (cleaned abstracts)
out_dir = paste(working_dir,"PubMed_PRUNE/",sep="")
result_dir = paste(working_dir,"PubMed_GENES/",sep="")
system(paste("mkdir ",result_dir,sep=""))

# all files in out_dir that have PRUNED in their name will be analyzed
file_list = list.files(out_dir)
files_ok = file_list[grep("PRUNED",file_list)]

### function for actual analysis of the abstracts (provided as a vector of text lines)
analyze_abstracts=function(raw.txt){
    # cycle through abstracts with gene-name extraction procedure
    for (j in 1:length(raw.txt)) { # takes a while... quite a big one
        txt=raw.txt[j]

        # these names have been retrieved by the translator
        replaced = unlist(strsplit(replacedf[j],","))
        replaced = replaced[which(replaced!="NA")]

        # remove all possible punctuation marks, but preserve the hyphen
        txt=paste(unlist(strsplit(txt,"-")),collapse="zzzzz") # "zubiquitilation of hyphens"
        txt=paste(unlist(strsplit(txt,"[[:punct:]]")),collapse=" ")
        txt=paste(unlist(strsplit(txt,"zzzzz")),collapse="-")

        # extract all remaining POTENTIAL gene names (assumptions are made below!)
        lst1=unlist(regmatches(txt, gregexpr('([[:punct:]]|[[:upper:]]|[0-9])+', txt)))
        lst2=unlist(regmatches(lst1, gregexpr('^[[:upper:]]{2,20}.*', lst1)))
        lst2 = sort(unique(lst2))
        lst2 = lst2[which(! lst2 %in% "REPLACED")]  # do not include artifacts

        # ELIMINATE GENES THAT ARE EASILY CONFUSED WITH TERMS FROM OBSTETRIC/MEDICAL/TECHNICAL FIELDS.
        # THE ACRONYMS LISTED MUST BE ADAPTED MANUALLY!
        lst3 = lst2[which(! lst2 %in% restricted_acronyms)]

        cumm1<<-c(cumm1,unique(lst3))
        cumm2<<-c(cumm2,unique(replaced))
        cumm3<<-c(cumm3,unique(c(lst3,replaced))) # will be used to estimate number of instances when a gene was mentioned
        rm(txt,lst1,lst2,lst3,replaced)
    }
}

### simple version
raw.txt=NULL
trfile="TRANSLATOR_misspelled_gene_names.txt"
translator = read.table(trfile, stringsAsFactors=F,h=T,sep="\t") # TRANSLATOR OF SOME (!) GENE NAMES
from_length = nchar(translator$from)
translator = translator[order(from_length,decreasing = T),]
write.table(translator,trfile,quote = F,row.names = F,col.names = T,sep="\t")

for(file_name in files_ok){
    print(paste("working on file",file_name))
    name_chunks = unlist(strsplit(file_name,"_"))
    outfile=paste(out_dir,"PubMed_webSearch_REPLACED",paste(name_chunks[4:8],collapse="_"),sep="_")
    outfile2=paste(out_dir,"PubMed_webSearch_REPLACED",paste(name_chunks[4:7],collapse="_"),"GENES",sep="_")

    ## look for genes found in translator table, using awk
    ## matching names are removed, so they wouldn't be counted twice
    print("translating gene names...")
    system(paste("awk -F'\t' 'FNR==NR && FNR>1{ a[$1]=$2 ; next }
                 {b=\"NA\"; n=0; for(i in a){
                            m=gsub(i,\"REPLACED\"); if(m>0){ b=b \",\" a[i] }
                        }
                   print b >\"",outfile2,"\"; print $0 >\"", outfile,"\"}' ",
                 trfile," ", paste(out_dir,file_name,sep=""),sep=""))

  raw.txt=readLines(outfile)
  replacedf=readLines(outfile2)
  print(paste("...translation complete. number of abstracts read:", length(raw.txt)))

  cumm1=NULL  # cummulation of potential gene names extracted from abstracts without TRANSLATOR
  cumm2=NULL  # cummulation of REAL gene names extracted from abstracts using TRANSLATOR
  cumm3=NULL  # cummulation of ALL unique-in-one-abstract gene names (ddtected with or without translator). for "times-mentioned" threshold

  # extract POTENTIAL genes (for example, all three-capital-letter words...)
  analyze_abstracts(raw.txt)

  # directly extracted POTENTIAL genes
  genes_1 = cumm1[which(cumm1 %in% hg_genes)]

  # genes detected via TRANSLATOR file
  genes_2 = cumm2[which(cumm2 %in% hg_genes)]

  # detected via BOTH METHODS
  genes_3 = cumm3[which(cumm3 %in% hg_genes)]

  ## attach counts and descriptions
  temp1 = data.frame(gene=names(table(genes_3)),freq=as.numeric(table(genes_3)),stringsAsFactors = F)
  temp2 = merge(temp1, hg, by.x="gene", by.y="HUGO", all.x=T); rm(temp1)
  gene.freq = temp2[rev(order(temp2$freq)),]; rm(temp2)

  ## write output
  result_file=paste("GENES",name_chunks[4],name_chunks[5],name_chunks[6],sep="_")
  write.table(gene.freq,paste(result_dir,result_file,sep=""),quote = F,row.names = F,col.names = F)
}
