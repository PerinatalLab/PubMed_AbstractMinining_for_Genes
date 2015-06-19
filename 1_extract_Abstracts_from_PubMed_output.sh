#!/bin/bash

# this script is controlled by a master-script written in R

in_file_name=$1
out_file_prfx=$2

# this directory is on WS2 computer (due to gawk)
data_dir=/home/jonas/Biostuff/MOBA_GESTAGE_GWAS/PREGNANCY_GENES/PubMed_2015Jun/

PubMed_results=${data_dir}/PubMed_RAW/${in_file_name} # downloaded from PubMed
PubMed_edited=${data_dir}/PubMed_DIGEST/temporary_file.txt
PubMed_digest_abstracts=${data_dir}PubMed_DIGEST/${out_file_prfx}_abstracts.txt
PubMed_digest_titles=${data_dir}PubMed_DIGEST/${out_file_prfx}_titles.txt


echo "editing references in the file..."
gawk 'BEGIN{RS="\n\n+"} {gsub("\n"," "); print $0}' ${PubMed_results} > ${PubMed_edited}

echo "writing out abstracts..."
gawk 'BEGIN{OFS="\t"; i=1} $1==i"."{pr=NR+3} NR==pr-1 && $0~/^\[Art/{pr++;next} \
NR==pr && $0~/^Author information/{pr++} NR==pr && $0~/^Collaborator/{pr++} \
NR==pr && $0~/^Erratum/{pr++} NR==pr && $0~/^Retraction /{i++; next} \
NR==pr && $0~/^Republished /{pr++} NR==pr && $0~/^Comment /{pr++} NR==pr{print i++, $0}' ${PubMed_edited} > ${PubMed_digest_abstracts}

echo "writing out titles..."
gawk 'BEGIN{OFS="\t"; i=1; t=""} $1==i"."{print (i++)-1, t; pr=NR+1} NR==pr{t=$0} $0~/^Retraction /{t=""}' ${PubMed_edited} | gawk '$2!=""'> ${PubMed_digest_titles}

rm ${PubMed_edited}

echo "...done!"
