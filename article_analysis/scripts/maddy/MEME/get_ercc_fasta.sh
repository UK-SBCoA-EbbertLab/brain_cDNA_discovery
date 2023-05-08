#!/bin/bash

singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/ERCC92/ERCC92.fa -bed ercc-00171.bed -s -name -fo ercc-00171.fa



