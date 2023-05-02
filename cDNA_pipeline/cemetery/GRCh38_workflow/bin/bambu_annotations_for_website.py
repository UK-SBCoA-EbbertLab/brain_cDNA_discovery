#!/usr/bin/env python

## Import libraries
import pandas as pd
import numpy as np
import sys
from datetime import date

TODAY = str(date.today())


'''
function name: parse_gtf_columns

purpose: parsing the last aggregate column of the gtf into useful columns and cleaning non-relevant columns

input: dataframe containining "raw" gtf

output: dataframe containing gtf with useful columns ["gene_id", "gene_name", "type"]
'''

def parse_gtf_columns(gtf, is_ref=True):

    ## Get gene ids
    gtf["gene_id"] = gtf['other'].str.split('"', expand=True)[1]

    ## Get transcript ids
    gtf["transcript_id"] = gtf["other"].str.split("transcript_id", expand=True)[1].str.split('"', expand=True)[1]

    if is_ref:

        ## Get gene names
        gtf["gene_name"] = gtf['other'].str.split('gene_name', expand=True)[1].str.split('"', expand=True)[1]

        ## Get gene types
        gtf["subtype"] = gtf["other"].str.split("gene_biotype", expand=True)[1].str.split('"', expand=True)[1]

        ## Get transcript names
        gtf["transcript_name"] = gtf["other"].str.split("transcript_name", expand=True)[1].str.split('"', expand=True)[1]

        ## Drop useless columns
        gtf.drop(columns=["dot_1", "dot_2", "other", "source", "strand", "start", "end", "type", "chr"], inplace=True)

    else:
        ## Drop useless columns
        gtf.drop(columns=["dot_1", "dot_2", "other", "source"], inplace=True)

    for col in gtf.columns:
        gtf.loc[gtf[col].isnull(), col] = np.NaN

    return gtf


'''
function name: merge_annotations

purpose: Merge useful/relevant information from both annotations while removing repeated and irrelevant information

input: Two different GTF annotations

output: One GTF annotation containing all the relevant information
'''

def merge_annotations(ref_gtf, bambu_gtf):

    ## Merge the two annotations
    merged_gtf = pd.merge(bambu_gtf, ref_gtf, on=['transcript_id'], how='left')
    merged_gtf["gene_id"] = merged_gtf["gene_id_x"]
    merged_gtf.drop(columns=["gene_id_x", "gene_id_y"], inplace=True)
    merged_gtf.drop_duplicates(inplace=True)

    ## Label novel transcripts
    merged_gtf.loc[merged_gtf["transcript_id"].str.startswith("tx."), "is_novel_transcript"] = True
    merged_gtf.loc[~merged_gtf["transcript_id"].str.startswith("tx."), "is_novel_transcript"] = False

    ## Label novel genes
    merged_gtf.loc[merged_gtf["gene_id"].str.startswith("gene."), "is_novel_gene"] = True
    merged_gtf.loc[~merged_gtf["gene_id"].str.startswith("gene."), "is_novel_gene"] = False

    ## Create temporary variable only containing novel transcripts
    temp = merged_gtf.loc[merged_gtf["is_novel_transcript"] == True]

    ## Annotate novel transcripts
    merged_tmp = pd.merge(temp, ref_gtf, on=['gene_id'], how='left')
    merged_tmp["transcript_id"] = merged_tmp["transcript_id_x"]
    merged_tmp["transcript_name"] = merged_tmp["transcript_name_x"]
    merged_tmp["subtype"] = merged_tmp["subtype_x"]
    merged_tmp["gene_name"] = merged_tmp["gene_name_y"]
    merged_tmp.drop(columns=["transcript_id_x", "transcript_id_y", "transcript_name_x",
                            "transcript_name_y", "subtype_x", "subtype_y", "gene_name_x", "gene_name_y"], inplace=True)
    merged_tmp.drop_duplicates(inplace=True)

    ## Return novel transcripts to original annotation
    merged_final = pd.merge(merged_gtf, merged_tmp, on=['chr', 'type', 'start', 'end', 'strand', 'transcript_id',
        'subtype', 'transcript_name', 'gene_id', 'is_novel_transcript', 'is_novel_gene'], how="left")

    ## Get gene names for novel transcripts of known genes
    merged_final.gene_name_x.fillna(merged_final.gene_name_y, inplace=True)
    merged_final["gene_name"] = merged_final["gene_name_x"]
    merged_final.drop(columns =["gene_name_x", "gene_name_y"], inplace=True)

    return merged_final


'''
function name: calculate_cpm

purpose: Calculate the CPM

input: File containing counts for every transcript in each sample

output: File containing CPM for every transcript in each sample
'''

def calculate_cpm(counts):

    ## Get the names of the count columns
    list_counts_cols = counts.columns[2:]

    new_col_names_list = ["transcript_id", "gene_id"]

    ## Loop through counts columns
    for col in list_counts_cols:
        
        ## Create new column name
        new_col_names_list.append(col.split("_nanopore")[0] + "_CPM")

        ## Get sum of counts column
        sum_col = counts[col].sum()

        ## Calculate CPM for new CPM column
        counts[col] = ((counts[col]/sum_col) * 1000000)

    ## Drop old Columns
    counts.columns = new_col_names_list

    return counts
 


def main():
    
    ## Load filenames from command lines
    ref_gtf_name = sys.argv[1]
    bambu_gtf_name = sys.argv[2]
    bambu_counts_name = sys.argv[3]
    output_name = sys.argv[4]



    ## Import reference GTF, only keep transcripts and exons, and parse through dataframe to get gene
    ## And transcript names
    ref_gtf = pd.read_csv(ref_gtf_name, delimiter="\t", header=4, low_memory=False,
                 names=["chr", "source", "type", "start", "end", "dot_1", "strand", "dot_2", "other"])
    
    ref_gtf = ref_gtf.loc[((ref_gtf["type"] == "transcript") | (ref_gtf["type"] == "exon"))]
    ref_gtf = parse_gtf_columns(ref_gtf, is_ref=True)
    
    ## Import Bambu GTF including novel genes/transcripts. Parse to data to extract gene_id and transcript_id
    bambu_gtf = pd.read_csv(bambu_gtf_name, header=None, delimiter="\t", low_memory=False,
                       names=["chr", "source", "type", "start", "end", "dot_1", "strand", "dot_2", "other"])
    
    bambu_gtf = bambu_gtf.loc[((bambu_gtf["type"] == "transcript") | (bambu_gtf["type"] == "exon"))]
    bambu_gtf = parse_gtf_columns(bambu_gtf, is_ref=False)
    
    ## Merge ref annotation and bambu annotation
    merged_gtf = merge_annotations(ref_gtf, bambu_gtf)
    
    ## Open bambu counts and calculate CPM
    bambu_counts = pd.read_csv(bambu_counts_name, delimiter="\t", low_memory=False, header=0)
    bambu_counts = calculate_cpm(bambu_counts)
    
    ## Create final transcripts annotation and save
    transcripts = merged_gtf.loc[merged_gtf["type"] == "transcript"]
    transcripts = pd.merge(transcripts, bambu_counts, on=["transcript_id", "gene_id"], how="left")
    transcripts_out_name = TODAY  + "_" + output_name + "_transcripts_annotation.csv"
    transcripts.to_csv(transcripts_out_name, index=False)
    
    ## Create final exons annotation and save
    exons = merged_gtf.loc[merged_gtf["type"] == "exon"]
    exons_out_name = TODAY  + "_" + output_name + "_exons_annotation.csv"
    exons.to_csv(exons_out_name, index=False)

main()

