#!/bin/bash

rclone copy --progress --copy-links --exclude="{*.fq,*.fastq}" --transfers 5 --checkers 5 \
    ./results/ \
    gemini32:/mnt/gemini3-4/PROJECTS/brain_discovery_cDNA/nextflow_output/
