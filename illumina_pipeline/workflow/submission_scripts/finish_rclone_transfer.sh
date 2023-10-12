#!/bin/bash

rclone copy --update --checksum --progress --copy-links --exclude="{*.fq,*.fastq}" --transfers 32 --checkers 32 \
        ./results/ gemini32:/mnt/gemini3-4/PROJECTS/brain_discovery_cDNA/nextflow_output/
