#!/bin/bash
#SBATCH --time=1-00:00:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=blat_novel_mito_alignments_%j    # Job name
#SBATCH --mem=100G                  ## Ask for 100G of RAM
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e blat_novel_mito_alignments_%j.err             # Error file for this job.
#SBATCH -o blat_novel_mito_alignments_%j.out             # Output file for this job.
#SBATCH -A coa_mteb223_uksr     # Project allocation account name (REQUIRED)


## BLAT novel genes against eachother
blat -minIdentity=80 -minScore=80 ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/novel_mt_transcripts.fa \
    ../../../../article_analysis/references/bernardo/Homo_sapiens.GRCh38_ERCC.fa new_spliced_mito_transcripts_alignments_to_GRCh38.psl

tail -n +6 new_spliced_mito_transcripts_alignments_to_GRCh38.psl > new_spliced_mito_transcripts_alignments_to_GRCh38_NO_HEADER.psl

./score_blat_output.awk new_spliced_mito_transcripts_alignments_to_GRCh38_NO_HEADER.psl > final_resutls_blat_alignments_spliced_mito_transcripts_to_GRCh38.txt


mv *.psl *GRCh38.txt ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/

