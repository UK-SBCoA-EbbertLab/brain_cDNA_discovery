#!/usr/bin/env python
# coding: utf-8


## Import libraries
import pandas as pd
import numpy as np
import sys
import os
import regex as re


## Load filenames from command lines
txt_path = sys.argv[1]
fastq_path = sys.argv[2]

## Create ouput file names
first_half_txt_out_name = "./" + re.sub("[a-zA-Z\d]*.txt", "firstHalf.txt", os.path.basename(txt_path))
first_half_fastq_out_name = "./" + re.sub("[a-zA-Z\d]*.fastq", "firstHalf.fastq", os.path.basename(fastq_path))

second_half_txt_out_name = "./" + re.sub("[a-zA-Z\d]*.txt", "secondHalf.txt", os.path.basename(txt_path))
second_half_fastq_out_name = "./" + re.sub("[a-zA-Z\d]*.fastq", "secondHalf.fastq", os.path.basename(fastq_path))

## Read fastq file line by line
fastq_in = open(fastq_path, 'r')
Lines = fastq_in.readlines()

## Create empty list for read ids in each file
list_read_id_first_half = []
list_read_id_second_half = []

## Open output fastq files for writing
fastq_first_half_out = open(first_half_fastq_out_name, 'w')
fastq_second_half_out = open(second_half_fastq_out_name, 'w')

length_lines = len(Lines)

## Create lists with read ids for split files and write lines to files
for i in range(0, length_lines, 8):

    list_read_id_first_half.append(Lines[i].split("@")[1].split(" ")[0])
    fastq_first_half_out.write(Lines[i])
    fastq_first_half_out.write(Lines[i+1])
    fastq_first_half_out.write(Lines[i+2])
    fastq_first_half_out.write(Lines[i+3])
    
    if ((i + 7) < length_lines):
        list_read_id_second_half.append(Lines[i+4].split("@")[1].split(" ")[0])
        fastq_second_half_out.write(Lines[i+4])
        fastq_second_half_out.write(Lines[i+5])
        fastq_second_half_out.write(Lines[i+6])
        fastq_second_half_out.write(Lines[i+7])
    
## Close files
fastq_first_half_out.close()
fastq_second_half_out.close()

## Set lines to None to save memory
Lines = None

## Read txt sequencing sammary file into dataframe
df_txt = pd.read_csv(txt_path, delimiter='\t', usecols=["read_id", "run_id", "channel", "start_time",
                 "sequence_length_template", "mean_qscore_template"], low_memory=False)

## Create txt files with corresponding reads to each of the fastq files created
df_txt_first_half = df_txt.loc[df_txt["read_id"].isin(list_read_id_first_half)]
df_txt_second_half = df_txt.loc[df_txt["read_id"].isin(list_read_id_second_half)]

## Write output sequencing summary files (pre and post wash)
df_txt_first_half.to_csv(first_half_txt_out_name, index=False, sep="\t")
df_txt_second_half.to_csv(second_half_txt_out_name, index=False, sep="\t")
