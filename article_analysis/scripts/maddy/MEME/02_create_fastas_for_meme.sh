#!/bin/bash

#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/five_prime_splice_sites.bed -name -fo ../../../data/maddy/MEME/five_prime_splice_sites.fa
#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/three_prime_splice_sites.bed -name -fo ../../../data/maddy/MEME/three_prime_splice_sites.fa
#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/exons.bed -name -fo ../../../data/maddy/MEME/exons.fa


#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/five_prime_splice_sites_nfk_12s.bed -s -name -fo ../../../data/maddy/MEME/five_prime_splice_sites_nfki_12s.fa
#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/three_prime_splice_sites_nfk_12s.bed -s -name -fo ../../../data/maddy/MEME/three_prime_splice_sites_nfk_12s.fa
#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/exons_nfk_12s.bed -s -name -fo ../../../data/maddy/MEME/exons_nfk_12s.fa


#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/five_prime_splice_sites_nfn_12s.bed -s -name -fo ../../../data/maddy/MEME/five_prime_splice_sites_nfn_12s.fa
#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/three_prime_splice_sites_nfn_12s.bed -s -name -fo ../../../data/maddy/MEME/three_prime_splice_sites_nfn_12s.fa
#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/exons_nfn_12s.bed -s -name -fo ../../../data/maddy/MEME/exons_nfn_12s.fa


singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../../references/bernardo/Homo_sapiens.GRCh38_ERCC.fa -bed ../../../data/maddy/MEME/five_prime_splice_sites_nfm_12s_filtered.bed -s -name -fo ../../../data/maddy/MEME/five_prime_splice_sites_nfm_12s_filtered.fa
singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../../references/bernardo/Homo_sapiens.GRCh38_ERCC.fa -bed ../../../data/maddy/MEME/three_prime_splice_sites_nfm_12s_filtered.bed -s -name -fo ../../../data/maddy/MEME/three_prime_splice_sites_nfm_12s_filtered.fa
singularity exec ../../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../../references/bernardo/Homo_sapiens.GRCh38_ERCC.fa -bed ../../../data/maddy/MEME/exons_nfm_12s_filtered.bed -s -name -fo ../../../data/maddy/MEME/exons_nfm_12s_filtered.fa

#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/five_prime_splice_sites_kfk_12s.bed -s -name -fo ../../../data/maddy/MEME/five_prime_splice_sites_kfk_12s.fa
#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/three_prime_splice_sites_kfk_12s.bed -s -name -fo ../../../data/maddy/MEME/three_prime_splice_sites_kfk_12s.fa
#singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_2022_09_28.sif bedtools getfasta -fi /project/mteb223_uksr/sequencing_resources/references/Ensembl/hg38_release_107/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed ../../../data/maddy/MEME/exons_kfk_12s.bed -s -name -fo ../../../data/maddy/MEME/exons_kfk_12s.fa



