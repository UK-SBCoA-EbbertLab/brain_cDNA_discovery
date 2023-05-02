#!/usr/bin/bash

fastq=$1
id=$2
cores=$3

isONclust  --t $cores  --ont --fastq $fastq \
             --outfolder clustering

isONclust write_fastq --N 1 --clusters clustering/final_clusters.tsv \
                      --fastq $fastq --outfolder clustering/fastq_files
