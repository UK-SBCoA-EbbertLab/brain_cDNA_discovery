#!/bin/bash

./OURS_GRCh38-107_discovery_subsample_0.2.sh

rclone copy --update --progress --copy-links --exclude="{*.fq,*.fastq}" --transfers 12 --checkers 12 \
    ./results/OURS_GRCh38-107_discovery_subsample_0.2/ \
    gemini32:/mnt/gemini3-4/PROJECTS/brain_discovery_cDNA/nextflow_output/OURS_GRCh38-107_discovery_subsample_0.2/

rclone copy --update --progress --copy-links --exclude="{*.fq,*.fastq}" --transfers 12 --checkers 12 \
    ./results/OURS_GRCh38-107_discovery_subsample_0.2/ \
    gemini32:/mnt/gemini3-4/PROJECTS/brain_discovery_cDNA/nextflow_output/OURS_GRCh38-107_discovery_subsample_0.2/

rclone copy --progress --update --copy-links --exclude="{*.fq,*.fastq}" --transfers 12 --checkers 12 \
    ./results/OURS_GRCh38-107_discovery_subsample_0.2/ \
    gemini32:/mnt/gemini3-4/PROJECTS/brain_discovery_cDNA/nextflow_output/OURS_GRCh38-107_discovery_subsample_0.2/

yes | rm -r ./results/ ./work/
