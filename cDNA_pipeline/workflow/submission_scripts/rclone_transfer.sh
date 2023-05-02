#!/bin/bash


for dir in ./results/*
do
    dir=$(basename "${dir}")

    rclone copy --transfers 64 --checkers 64 --copy-links --progress \
        "./results/${dir}/" "gemini1-2:/mnt/gemini1-6/mteb223_uksr/PROCESSED_DATA/bag222/brain_discovery/data/raw/${dir}/"

done
