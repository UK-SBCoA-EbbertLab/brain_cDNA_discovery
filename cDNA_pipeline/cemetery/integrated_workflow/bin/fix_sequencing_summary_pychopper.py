#!/usr/bin/env python

## Import libraries
import pandas as pd
import sys

## Load filenames from command lines
fastq_filename = sys.argv[1]
txt_filename = sys.argv[2]
output_name = sys.argv[3]

## Read fastq file line by line
fastq = open(fastq_filename, 'r')
Lines = fastq.readlines()

fastq_read_ids = []

# Strips the newline character
for line in Lines:

    ## Get line content
    line_content = line.strip()

        ## Get read_id from fastq_file
    if ((line_content[0] == "@") & (line_content[-8:-2] == "strand")):
        read_id = line_content.split("|")[1].split(' ')[0]
        fastq_read_ids.append(read_id)

## Read txt sequencing sammary file in
df_txt = pd.read_csv(txt_filename, delimiter='\t')

## Remove useless columns
df_txt = df_txt.drop(df_txt.columns.difference(["read_id", "run_id",
                                                "channel", "start_time", "sequence_length_template",
                                                "mean_qscore_template"]), axis=1, inplace=False)

## Remove points where the txt files were joined (Headers in the middle of the data)
df_txt = df_txt.loc[df_txt["channel"]!="channel"]
df_txt[["channel", "start_time", "sequence_length_template",
        "mean_qscore_template"]] = df_txt[["channel", "start_time", "sequence_length_template",
                                       "mean_qscore_template"]].apply(pd.to_numeric, errors="ignore")
df_txt.dropna(inplace=True)

## Only keep reads in the sequencing summary that were kept after pychopper filtering
df_txt = df_txt.loc[df_txt["read_id"].isin(fastq_read_ids)]

## Save new filtered sequencing summary file as a txt file
df_txt.to_csv(output_name, index=False, sep="\t")
