#!/bin/bash

## scripts to convert and annotate sites using genomation 

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -c|--chromolist)
            chromolist="$2"
            shift # past argument
            shift # past value
            ;;
        -o|--organism)
            organism="$2"
            shift # past argument
            shift # past value
            ;;
        -s|--sites)
            sites="$2"
            shift # past argument
            shift # past value
            ;;
        -g|--genome)
            genome="$2"
            shift # past argument
            shift # past value
            ;;
        *)  # unknown option
            echo "Unknown option: $1"
            shift # past argument
            ;;
    esac
done

# Check if any required arguments are missing
if [ -z "$chromolist" ] || [ -z "$organism" ] || [ -z "$sites" ] || [ -z "$genome" ]; then
    echo "Usage: $0 -c|--chromolist <chromosome_list_file> -o|--organism <organism_name> -s|--sites <sites_file> -g|--genome <genome_file>"
    exit 1
fi

#### get chromosomes as number

cat ${chromolist} | while read line
do
  code=$(echo $line | cut -d " " -f1)
  number=$(echo $line | cut -d " " -f2)
  echo "s/${code}/${number}/"
done  > chr.sed.script

### converting csv to bed ####

cat ${sites} | sed 's/,/\t/g' | awk 'NR > 1 {print $2,$3,$3}' | sed 's/ * /\t/g' | sed 's/"//g' > ${organism}.sites.bed

### converting gtf files to bed12 and right chromosome codes ####

./gtfToGenePred -ignoreGroupsWithoutExons ${genome} genes.pred
./genePredToBed genes.pred ${organism}.ann.bed

rm genes.pred

sed -f chr.sed.script ${organism}.ann.bed > ${organism}.coded.ann.bed

rm ${organism}.ann.bed

echo "file conversions done"

#### R part to annotate with genomation (building R script) ####

cat > genomation_run.r << EOF

## loading libraries

if (!requireNamespace("genomationData", quietly = TRUE)) install.packages("genomationData", repos = "https://cloud.r-project.org/")
if (!requireNamespace("genomation", quietly = TRUE)) install.packages("genomation", repos = "https://cloud.r-project.org/")
if (!requireNamespace("GenomicRanges", quietly = TRUE)) install.packages("GenomicRanges", repos = "https://cloud.r-project.org/")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = "https://cloud.r-project.org/")

library(genomationData)
library(GenomicRanges)
library(readr)
library(genomation)

## loading data

## cpg sites

#coord = read.table("${organism}.sites.bed")
#cpgsites = readBed("${organism}.sites.bed")

coord = read.csv("${sites}")

cpgsites = readBed("${organism}.sites.bed")

## annotated genome

genomeann = readTranscriptFeatures("${organism}.coded.ann.bed")
 
## running annotation

annot.list = annotateWithGeneParts(cpgsites, genomeann)

## get distance from TSS and gene ids

tssdist = getAssociationWithTSS(annot.list)

position = getMembers(annot.list)

output <- cbind(coord,tssdist,position)

tooutput <- output[, c(1, 2, 5, 6, 7, 8, 9,10)]

write.table(tooutput, "${organism}_annotated_sites.txt", sep = "\t", row.names = FALSE, quote = FALSE)

EOF

module load R

Rscript genomation_run.r

echo "annotation done"

### getting transcript to gene information 

cat ${genome} | cut -d "\"" -f2,4 | sed 's/\"/\t/g' | awk -v OFS="\t" '!/^#/ && NF == 2 {print $2,$1}' | sort | uniq  > ${organism}_gene_transcript_info


## running association of transcripts to gene 

cat > gene_names.py << EOF

import pandas as pd

existing_table = pd.read_table('${organism}_annotated_sites.txt')
info_to_append = pd.read_table('${organism}_gene_transcript_info', header=None, names=['Column1', 'Column2'])
merged_table = pd.merge(existing_table, info_to_append, left_on='feature.name', right_on='Column1', how='left')
merged_table.rename(columns={'Column2': 'gene.name'}, inplace=True)
merged_table.rename(columns={'V1': 'chromossome'}, inplace=True)
merged_table.rename(columns={'V2': 'site'}, inplace=True)
merged_table.drop(columns=['Column1'], inplace=True)

merged_table.to_csv('${organism}_annotated_sites_complete.txt', sep='\t', index=False)

EOF

echo "running final python steps"

python gene_names.py 

echo "script done final annotation file : ${organism}_annotated_sites_complete.txt" 