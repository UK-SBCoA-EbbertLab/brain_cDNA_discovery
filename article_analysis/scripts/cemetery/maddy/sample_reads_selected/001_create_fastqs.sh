#!/bin/bash

samples=$1

echo $samples

for bam in ${samples[@]}
do
	file_name="$(basename $bam)"
	samtools view -h $bam 'GL000214.1:125300-133400'  | samtools fastq > "/pscratch/mteb223_uksr/brain_cDNA_discovery/article_analysis/data/maddy/sample_reads_selected/${file_name}_GL000214_125300_133400_region.fastq"
	samtools view -h $bam 'GL000214.1'  | samtools fastq > "/pscratch/mteb223_uksr/brain_cDNA_discovery/article_analysis/data/maddy/sample_reads_selected/${file_name}_GL000214_region.fastq"
done


