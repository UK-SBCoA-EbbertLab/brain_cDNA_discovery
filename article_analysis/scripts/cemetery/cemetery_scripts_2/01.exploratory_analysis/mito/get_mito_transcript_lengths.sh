#!/bin/bash

awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/novel_mt_transcripts.fa \
    > ../../../data/bernardo/processed/05.mitochondrial_novel_transcripts_probing/novel_mt_transcript_lengths.txt
