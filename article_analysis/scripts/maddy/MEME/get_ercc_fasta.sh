#!/bin/bash

singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi \
    /project/mteb223_uksr/sequencing_resources/references/ERCC92/ERCC92.fa -bed ercc-00171.bed -s -name -fo ercc-00171.fa



