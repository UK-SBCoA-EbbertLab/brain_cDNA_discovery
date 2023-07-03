#!/bin/bash


sbatch mafft_alignment.sbatch ../../../data/bernardo/processed/05.sanger_sequencing/sanger_sequencing_fasta_files/MAFFT_alignment_BambuTx1845_forward.fa \
    ../../../data/bernardo/processed/05.sanger_sequencing/mafft_alignment_results/MAFFT_alignment_BambuTx1845_forward_OUTPUT.txt

sbatch mafft_alignment.sbatch ../../../data/bernardo/processed/05.sanger_sequencing/sanger_sequencing_fasta_files/MAFFT_alignment_BambuTx1845_reverse.fa \
    ../../../data/bernardo/processed/05.sanger_sequencing/mafft_alignment_results/MAFFT_alignment_BambuTx1845_reverse_OUTPUT.txt

sbatch mafft_alignment.sbatch ../../../data/bernardo/processed/05.sanger_sequencing/sanger_sequencing_fasta_files/MAFFT_alignment_BambuTx2703_forward.fa \
    ../../../data/bernardo/processed/05.sanger_sequencing/mafft_alignment_results/MAFFT_alignment_BambuTx2703_forward_OUTPUT.txt

sbatch mafft_alignment.sbatch ../../../data/bernardo/processed/05.sanger_sequencing/sanger_sequencing_fasta_files/MAFFT_alignment_BambuTx2703_reverse.fa \
    ../../../data/bernardo/processed/05.sanger_sequencing/mafft_alignment_results/MAFFT_alignment_BambuTx2703_reverse.fa
