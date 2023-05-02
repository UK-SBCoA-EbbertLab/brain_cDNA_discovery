#!/bin/bash


#singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif ./AD_transcript_boxplots_AD_vs_CT.R
#singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif ./AD_transcript_boxplots_M_vs_F.R
singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif ./results_boxplots_AD_gene.R
singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif ./results_boxplots_AD_transcript.R
singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif ./results_boxplots_SEX_gene.R
singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif ./results_boxplots_SEX_transcript.R

