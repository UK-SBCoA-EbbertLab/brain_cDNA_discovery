#!/bin/bash

#sbatch --wait minimap_index.sh

files="/pscratch/mteb223_uksr/brain_cDNA_discovery/article_analysis/data/maddy/sample_reads_selected/*.fastq"

for fastq in ${files[@]}
do
	sbatch align_to_chm13.sh chm13_ref.mmi $fastq

done
