#!/bin/bash

sleep 6h

rclone copy --progress --copy-links --exclude="{*.fq,*.fastq}" --transfers 12 --checkers 12 \
    ./results/OURS_GRCh38-88_discovery/ \
    gemini32:/mnt/gemini3-4/PROJECTS/brain_discovery_cDNA/nextflow_output/OURS_GRCh38-88_discovery/

rclone sync --progress --copy-links --exclude="{*.fq,*.fastq}" --transfers 12 --checkers 12 \
    ./results/OURS_GRCh38-88_discovery/ \
    gemini32:/mnt/gemini3-4/PROJECTS/brain_discovery_cDNA/nextflow_output/OURS_GRCh38-88_discovery/

rclone sync --progress --copy-links --exclude="{*.fq,*.fastq}" --transfers 12 --checkers 12 \
    ./results/OURS_GRCh38-88_discovery/ \
    gemini32:/mnt/gemini3-4/PROJECTS/brain_discovery_cDNA/nextflow_output/OURS_GRCh38-88_discovery/


yes | rm -r ./results/ ./work/
