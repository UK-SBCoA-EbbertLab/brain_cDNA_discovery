#!/bin/bash


## New from Known
singularity exec ../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../references/Homo_sapiens.GRCh38_ERCC.fa -bed ../../data/processed/MEME/five_prime_splice_sites_nfk_12s.bed -s -name -fo ../../data/processed/MEME/five_prime_splice_sites_nfk_12s.fa
singularity exec ../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../references/Homo_sapiens.GRCh38_ERCC.fa -bed ../../data/processed/MEME/three_prime_splice_sites_nfk_12s.bed -s -name -fo ../../data/processed/MEME/three_prime_splice_sites_nfk_12s.fa


## New from New
singularity exec ../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../references/Homo_sapiens.GRCh38_ERCC.fa -bed ../../data/processed/MEME/five_prime_splice_sites_nfn_12s.bed -s -name -fo ../../data/processed//MEME/five_prime_splice_sites_nfn_12s.fa
singularity exec ../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../references/Homo_sapiens.GRCh38_ERCC.fa -bed ../../data/processed/MEME/three_prime_splice_sites_nfn_12s.bed -s -name -fo ../../data/processed/MEME/three_prime_splice_sites_nfn_12s.fa


## New from Mitochondrial
singularity exec ../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../references/Homo_sapiens.GRCh38_ERCC.fa -bed ../../data/processed/MEME/five_prime_splice_sites_nfm_12s.bed -s -name -fo ../../data/processed//MEME/five_prime_splice_sites_nfm_12s.fa
singularity exec ../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../references/Homo_sapiens.GRCh38_ERCC.fa -bed ../../data/processed/MEME/three_prime_splice_sites_nfm_12s.bed -s -name -fo ../../data/processed/MEME/three_prime_splice_sites_nfm_12s.fa


## Known from Known
singularity exec ../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../references/Homo_sapiens.GRCh38_ERCC.fa -bed ../../data/processed/MEME/five_prime_splice_sites_kfk_12s.bed -s -name -fo ../../data/processed//MEME/five_prime_splice_sites_kfk_12s.fa
singularity exec ../../../singularity_containers/bernardo_article_analysis.sif bedtools getfasta -fi ../../references/Homo_sapiens.GRCh38_ERCC.fa -bed ../../data/processed/MEME/three_prime_splice_sites_kfk_12s.bed -s -name -fo ../../data/processed/MEME/three_prime_splice_sites_kfk_12s.fa


