#!/bin/bash


rclone copy --progress --copy-links --transfers 64 --checkers 64 ./results/ gemini32:/mnt/gemini3-4/PROJECTS/brain_discovery_cDNA/
