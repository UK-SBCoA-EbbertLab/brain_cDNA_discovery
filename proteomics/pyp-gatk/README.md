# Protein Fasta Creation Scripts

The scripts in this directory were used to create the protein fasta file used in the fragpipe mass spec analysis. 

### Scripts

1. `create_uniprot_peptide_fasta.sh` -- This is the main script to run that will call the other scripts. It will grab the transcript ids and gene ids from the NEW and KNOWN gtfs, grab the protein sequence from biomaRt for the known transcripts, and generate protein sequence (in 3 potential reading frames) for the new transcripts using pyp-gatk.
2. `create_biomaRt_annotated_database.R` -- This script is passed the transcript ids for the known transcripts from genes where a new transcript was identified and uses biomaRt to get the appropriate protein sequence. 
3. `create_uniprot_headers.py` -- This script is used to create a uniprot style header for the fasta file.

