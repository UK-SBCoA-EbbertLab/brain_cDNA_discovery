#!/bin/bash

singularity exec ../../../singularity_containers/bernardo_article_analysis.sif ./01_create_new_gtf_for_fastas.R

./02_create_fastas_for_meme.sh



