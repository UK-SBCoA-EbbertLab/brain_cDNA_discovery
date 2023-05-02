#!/usr/bin/bash

fastq=$1
id=$2
cores=$3
i=$4

run_isoncorrect --fastq_folder $fastq  --outfolder correction/ --split_mod 40 --residual $i --t $cores

touch "./${id}_${i}_corrected_reads.fq"

for f in 'correction/*/corrected_reads.fastq';
do
    echo $f
    cat $f >> "${id}_${i}_corrected_reads.fq"
done
