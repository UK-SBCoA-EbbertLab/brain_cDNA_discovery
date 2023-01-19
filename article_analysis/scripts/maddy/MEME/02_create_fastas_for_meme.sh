#!/bin/bash
singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/five_prime_splice_sites.bed -name -fo ../../../data/maddy/MEME/five_prime_splice_sites.fa

singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/three_prime_splice_sites.bed -name -fo ../../../data/maddy/MEME/three_prime_splice_sites.fa

singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/exons.bed -name -fo ../../../data/maddy/MEME/exons.fa
