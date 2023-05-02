#!/usr/bin/env python

## Import libraries
import pandas as pd
import sys

## Load filenames from command lines
txt_filename = sys.argv[1]
output_name = sys.argv[2]

## Read txt sequencing sammary file in
df_txt = pd.read_csv(txt_filename, delimiter='\t')

## Remove useless columns
df_txt = df_txt[["read_id", "run_id", "channel", "start_time", "sequence_length_template", "mean_qscore_template"]].copy()

## Remove points where the txt files were joined (Headers in the middle of the data)
df_txt = df_txt.loc[df_txt["channel"]!="channel"]
df_txt[["channel", "start_time", "sequence_length_template",
        "mean_qscore_template"]] = df_txt[["channel", "start_time", "sequence_length_template",
                                       "mean_qscore_template"]].apply(pd.to_numeric, errors="ignore")
df_txt.dropna(inplace=True, axis=0)

## Save new filtered sequencing summary file as a txt file
df_txt.to_csv(output_name, index=False, sep="\t")
