#!/usr/bin/env python

## Import libraries
import numpy as np
import pandas as pd
import sys

## Load filenames from command lines
fastq_filename = sys.argv[1]
txt_filename = sys.argv[2]
output_name = sys.argv[3]

## Read fastq file line by line
fastq = open(fastq_filename, 'r')
Lines = fastq.readlines()

original_read_ids_list = []
pychopper_read_ids_list = []
read_length_list = []

# Strips the newline character
for line in Lines:

    ## Get line content
    line_content = line.strip()

    ## Get read_ids from fastq_file
    if ((line_content[0] == "@") & (line_content.find("runid=") != -1)):
        
        ## Get original and pychopper read ids
        pychopper_read_id = line_content.split('@')[1].split(' ')[0]
        original_read_id = line_content.split("|")[1].split(" ")[0]
        
        ## Append ids to their respective lists
        original_read_ids_list.append(original_read_id)
        pychopper_read_ids_list.append(pychopper_read_id)
        is_read_id = True

    
    ## If this is a sequence line
    elif ((line_content[0] in ["C", "T", "G", "A"]) & (line_content.isalpha()) & (is_read_id)):
        
        ## Get new read size after pychopper
        length = len(line_content)
        read_length_list.append(length)
        is_read_id = False


## Read txt sequencing sammary file in
df_txt = pd.read_csv(txt_filename, delimiter='\t')

## Remove useless columns
df_txt = df_txt[["read_id", "run_id", "channel", "start_time", "sequence_length_template", "mean_qscore_template"]].copy()

## Remove points where the txt files were joined (Headers in the middle of the data)
df_txt = df_txt.loc[df_txt["channel"]!="channel"].copy()
df_txt[["channel", "start_time", "sequence_length_template",
        "mean_qscore_template"]] = df_txt[["channel", "start_time", "sequence_length_template",
                                       "mean_qscore_template"]].apply(pd.to_numeric, errors="ignore").copy()
df_txt.dropna(inplace=True, axis=0)

## Create dataframe with pychopper info
df_read_converter = pd.DataFrame()
df_read_converter["read_id"] = np.asarray(original_read_ids_list)
df_read_converter["pychopper_read_id"] = np.asarray(pychopper_read_ids_list)
df_read_converter["pychopper_sequence_length_template"] = np.asarray(read_length_list).astype('int32')


## Merge with original txt and substitute old columns
df_txt_final = df_txt.merge(df_read_converter, on="read_id", how="outer")
df_txt_final["pychopper_read_id"].fillna(df_txt_final["read_id"], inplace=True)
df_txt_final["pychopper_sequence_length_template"].fillna(df_txt_final["sequence_length_template"], inplace=True)
df_txt_final["read_id"] = df_txt_final["pychopper_read_id"].copy()
df_txt_final["sequence_length_template"] = df_txt_final["pychopper_sequence_length_template"].copy()
df_txt_final.drop(columns=["pychopper_read_id", "pychopper_sequence_length_template"], inplace=True)


## Save new filtered sequencing summary file as a txt file
df_txt_final.to_csv(output_name, index=False, sep="\t")
