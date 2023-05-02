#!/usr/bin/env python

## Import libraries
import pandas as pd
import sys
import csv
import numpy as np


'''
function name: parse_df_columns

purpose: parsing the last aggregate column of the gtf/gff3 into useful columns and cleaning non-relevant columns

input: dataframe containining "raw" gtf/gff3 data

output: dataframe containing gtf with useful columns
'''

def parse_df_columns(df, is_ref=True):
    
    if is_ref:
        
        ## Get gene ids
        df["gene_id"] = df["other"].str.split("source_gene=", expand=True)[1].str.split(';', expand=True)[0]

        ## Get transcript ids
        df["transcript_id"] = df["other"].str.split("source_transcript=", expand=True)[1].str.split(';', expand=True)[0]
        
        ## Get CHM gene_ids 
        df["CHM_gene_id"] = df["other"].str.split("gene_id=", expand=True)[1].str.split(';', expand=True)[0]

        ## Get transcript ids
        df["CHM_transcript_id"] = df["other"].str.split("transcript_id=", expand=True)[1].str.split(';', expand=True)[0]
        
        ## Only keep relevant 
        df.drop(columns="other", inplace=True)
        
        ## Drop duplicates
        df.drop_duplicates(inplace=True)
        
    for col in df.columns:
        df.loc[df[col].isnull(), col] = np.NaN
    
    return df


def main():

    ## Define file names
    gff_name = sys.argv[1]
    output_name = sys.argv[2]

    ## Open gff reference file
    gff = pd.read_csv(gff_name, delimiter="\t", header=1,
                                 names = ["chr", "source", "type", "start", "end", "dot1", "strand", "dot2", "other"])
    ## Only keep transcripts
    gff = gff.loc[gff["type"].isin(["transcript", "exon"])]

    ## Parse through "other" column to extract important information
    gff = parse_df_columns(gff, is_ref=True)

    ## Change name of duplicate Ensembl IDs to CHM IDs
    gff.loc[gff["transcript_id"] == "N/A", "transcript_id"] = gff["CHM_transcript_id"]
    gff_transcripts = gff.loc[gff["type"] == "transcript"].copy()
    gff_transcripts = gff_transcripts[["transcript_id", "CHM_transcript_id"]].drop_duplicates()
    gff_transcripts = gff_transcripts[gff_transcripts['transcript_id'].duplicated() == True]
    dup_trans = gff_transcripts["transcript_id"].dropna().values.tolist()
    gff.loc[gff["transcript_id"].isin(dup_trans), "transcript_id"] = gff["transcript_id"] + "(" + gff["CHM_transcript_id"] + ")"

    ## Change name of duplicate gene ids to CHM ids
    gff.loc[gff["gene_id"] == "None", "gene_id"] = gff["CHM_gene_id"]
    gff_genes = gff.loc[gff["type"] == "transcript"].copy()
    gff_genes = gff_genes[["gene_id", "CHM_gene_id"]].drop_duplicates()
    gff_genes = gff_genes[gff_genes['gene_id'].duplicated() == True]
    dup_genes = gff_genes["gene_id"].dropna().values.tolist()
    gff.loc[gff["gene_id"].isin(dup_genes), "gene_id"] = gff["gene_id"] + "(" + gff["CHM_gene_id"] + ")"


    ## Create other column with proper transcript ids and gene ids
    gff["other"] = 'gene_id "' + gff["gene_id"] + '"; ' \
                        + 'transcript_id "' + gff["transcript_id"] + '";'

    ## Drop columns that don't belong in gtf file
    gff.drop(columns = ["gene_id", "transcript_id", "CHM_gene_id", "CHM_transcript_id"], inplace=True)

    ## Sort by index
    gff.sort_index(inplace=True)

    ## Save file
    gff.to_csv(output_name, index=False, header=False, sep="\t", quoting=csv.QUOTE_NONE)
        
main()
