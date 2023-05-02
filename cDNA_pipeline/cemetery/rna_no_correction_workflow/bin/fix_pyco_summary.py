#!/usr/bin/env python

import pandas as pd
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, delimiter='\t')


df = df.drop(df.columns.difference(["read_id", "run_id", "channel", "start_time", "sequence_length_template", "mean_qscore_template"]), axis=1, inplace=False)

df = df.loc[df["channel"]!="channel"]
df[["channel", "start_time", "sequence_length_template", "mean_qscore_template"]] = df[["channel", "start_time", "sequence_length_template", "mean_qscore_template"]].apply(pd.to_numeric, errors="ignore")
df.dropna(inplace=True)



df.to_csv(output_file, index=False, sep="\t")
