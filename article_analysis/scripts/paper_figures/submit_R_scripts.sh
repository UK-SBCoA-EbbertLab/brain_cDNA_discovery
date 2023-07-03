#!/bin/bash


singularity exec ../../../singularity_containers/bernardo_article_analysis.sif ./make_volcano_plots_AD_vs_CT_no_intron.R

singularity exec ../../../singularity_containers/bernardo_article_analysis.sif ./results_boxplots_AD_gene_no_intron.R

singularity exec ../../../singularity_containers/bernardo_article_analysis.sif ./results_boxplots_AD_transcript_no_intron.R
