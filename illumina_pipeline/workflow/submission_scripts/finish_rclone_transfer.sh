#!/bin/bash

for i in {1..10}
do

   rclone copy --update --progress --copy-links --exclude="{*.fq,*.fastq}" --transfers 32 --checkers 32 \
           ./results/ROSMAP_illumina_DorsoLateralPreFrontalCortex_UNIQUE/ \
               gemini32:/mnt/gemini3-4/PROJECTS/brain_discovery_cDNA/nextflow_output/ROSMAP_illumina_DorsoLateralPreFrontalCortex_UNIQUE/

done
