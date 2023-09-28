#!/usr/bin/env python
# coding: utf-8


## Import libraries
import pandas as pd
import numpy as np
import sys
import os
import regex as re
import random
from io import StringIO


"""
Subsample a FASTQ file represented as a list of strings by a given percentage.
    
Parameters:
- file_list (list of strings): The FASTQ file content represented as a list of strings.
- percentage (float): The desired subsampling percentage (0-100).
    
Returns:
- fastq_list (list of strings): The subsampled FASTQ content as a list of strings.
- read_ids (list): A list of read IDs in the subsampled FASTQ content (without '@' symbol).
"""

def subsample_fastq(file_list, percentage):
    
    # Group the file_list into individual reads (assuming each read consists of 4 lines)
    entries = [file_list[i:i+4] for i in range(0, len(file_list), 4)]

    # Randomly select the desired number of reads
    n = int(len(entries) * (percentage / 100))
    sampled = random.sample(entries, n)

    # Extract read IDs and create the subsampled FASTQ content list
    read_ids = []
    fastq_list = []
    for entry in sampled:
        read_id = entry[0][1:].split(" ")[0]  # Remove '@' from the ID
        read_ids.append(read_id)
        fastq_list.extend(entry)

    return fastq_list, read_ids


def main():

    ## Load filenames from command lines
    txt_path = sys.argv[1]
    fastq_path = sys.argv[2]

    ## Read fastq file line by line
    fastq_in = open(fastq_path, 'r')
    Lines = fastq_in.readlines()

    ## Read txt sequencing sammary file into dataframe
    df_txt = pd.read_csv(txt_path, delimiter='\t', usecols=["read_id", "run_id", "channel", "start_time",
                    "sequence_length_template", "mean_qscore_template"], low_memory=False)

    ## Loop over subsampling range
    for subsample in range(5, 101, 5):

        ## Define subsample percent
        subsample_percent = (subsample/100)

        ## Create downsample string
        down_txt = "_downsampled-" + str(subsample/100) + ".txt"
        down_fastq = "_downsampled-" + str(subsample/100) + ".fastq"

        ## Create ouput file names
        txt_out_name = "./" + re.sub(".txt", down_txt, os.path.basename(txt_path))
        fastq_out_name = "./" + re.sub(".fastq", down_fastq, os.path.basename(fastq_path))

        ## Open output fastq files for writing
        fastq_out = open(fastq_out_name, 'w')

        ## Create subsampled output and list of read_ids in it
        subsampled_lines, list_read_ids = subsample_fastq(Lines, subsample) 

        print(list_read_ids)

        ## Write the randomly subsampled output to a new fastq file
        fastq_out.writelines(subsampled_lines)
        
        ## Close fastqfile
        fastq_out.close()

        ## Create txt files with corresponding reads to the fastq file created
        df_txt_out = df_txt.loc[df_txt["read_id"].isin(list_read_ids)].copy()

        ## Write output sequencing summary file
        df_txt_out.to_csv(txt_out_name, index=False, sep="\t")

if __name__ == "__main__":
    main()
