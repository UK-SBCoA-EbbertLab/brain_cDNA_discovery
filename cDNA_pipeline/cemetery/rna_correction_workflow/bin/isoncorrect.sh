#!/usr/bin/bash

fastq=$1
id=$2
cores=$3


run_isoncorrect --t $cores --fastq_folder $fastq  --outfolder ./correction/


touch "./${id}_corrected_reads.fq"

for f in 'correction/*/corrected_reads.fastq';
do
    echo $f
    cat $f >> "${id}_all_corrected_reads.fq"
done
