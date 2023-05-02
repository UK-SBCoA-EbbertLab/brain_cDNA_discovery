#!/usr/bin/bash

fastq=$1
id=$2
cores=$3

isONclust  --t $cores  --ont --fastq $fastq \
             --outfolder clustering

isONclust write_fastq --N 1 --clusters clustering/final_clusters.tsv \
                      --fastq $fastq --outfolder clustering/fastq_files


run_isoncorrect --t $cores --fastq_folder clustering/fastq_files  --outfolder ./correction/


touch "./${id}_corrected_reads.fq"

for f in 'correction/*/corrected_reads.fastq';
do
    echo $f
    cat $f >> "${id}_all_corrected_reads.fq"
done
