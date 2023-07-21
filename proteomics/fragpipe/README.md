# fragpipe

## Getting Started

### 1) Install Fragpipe and necessary dependencies[^1].

 - Information on how Fragpipe works and how to install can be found here [here](https://fragpipe.nesvilab.org/).

### 2) Download public proteomics (mass spec) data.

 - Proteomics (Mass spec) data from cell-lines used in this experiment are publicly available [here](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024364). For more information about this data see: https://pubmed.ncbi.nlm.nih.gov/36959352/
	
 - Proteomics (Mass spec) data from round 2 of the ROSMAP TMT brain Proteomics are puclicly available [here](https://www.synapse.org/#!Synapse:syn17015098). For more information about this data see: https://www.nature.com/articles/s41597-020-00650-8

### 3) Clone this repository using the command below:

```
	git clone https://github.com/UK-SBCoA-EbbertLab/brain_cDNA_discovery
```

### 4) Modify `fragpipe.workflow`, and `fragpipe_submission_script`. The files used for the original submission can be found on `./brain_tissue_submission/` and `./cell_line_submission/`

 - Use X11 Forwarding to load the workflow into the Fragpipe GUI and make necessary changes.
 
 - Modify the submission script to include the correct paths and working directory. 

 - If using a workload manager then give the job enough time to run (36 hours should be more than enough).

## Downstream Analysis

### Cell-lines

	1) Open the `combined_peptide.tsv` in Excel.
	
	2) Filter the "Protein" column for IDs containing "Bambu" and filter the "Mapped Genes" and "Mapped Proteins" columns for blank entries.

	3) Run `../scripts/count_peptides.sh` with the remaining peptide sequences and `peptide_database_uniprot_header_run_MEDIAN_cpm_greater_than_one.fa` to get the number of times each peptide sequence shows up in the database.

	4) Filter out any peptide sequences that had a count > 1.
	
	5) Create a column to sum each spectral count entry for each peptide from each sample.

	6) Group by "Protein" and sum total number of spectral count for proteins.

	6) Filter Protein level total spectral count for values > 5.

	8) Use the remaining protein IDs to get corresponding list of peptide sequences.

	9) Run `../scripts/create_fasta.sh` with the list of protein IDs and peptide sequences to produce a fasta file.

 	10) Perform a blast protein search using the fasta file generated with `create_fasta.sh`. The parameters we used for our blast search can be found on [Zenodo](XXX).

  	11) Remove any peptides that have a blast hit with 100% identity match and 100% query coverage from the file created in step 5.

   	12) Run steps 6-9 again for the newly created file from step 11.

	13) These proteins with > 5 unique hits are considered validated in the analysis done in the article.

### Brain
	
	Since there are 14 batches for this data, downstream analysis is slightly different.

	1) Modify `../scripts/Merge_Peptides.ipynb` with the correct directory for the peptide.tsv files.

	2) Open `../scripts/Merge_Peptides.ipynb` as a jupyter notebook, run,  and export the column of peptide sequences as a csv file.

	3) Run `../scripts/count_peptides` with the peptide sequences csv file that was just created and `peptide_database_uniprot_header_run_MEDIAN_cpm_greater_than_one.fa` to get the number or times each peptide sequence shows up in the database.

	4) Filter out any peptide sequence that had a count > 1.

 	5) Group by "Protein" and sum total number of spectral count for proteins.

	6) Filter Protein level "Spectral Count" for values > 5.
	
	7) Run `../scripts/create_fasta.sh` with  the remaining protein IDs and corresponding peptide sequences to produce a fasta file.

 	8) Perform a blast protein search using the fasta file generated with `create_fasta.sh`. The parameters we used for our blast search can be found on [Zenodo](XXX).

  	9) Remove any peptides that have a blast hit with 100% identity match and 100% query coverage from the file created in step 5.

   	10) Run steps 5 and 6 again for the newly created file from step 8.

	11) These proteins with > 5 unique hits are considered validated in the analysis done in the article.

[^1]:This will include Java 9+, [MSFragger](https://msfragger.nesvilab.org/), [Philosopher](https://philosopher.nesvilab.org/), [IonQuant](http://ionquant.nesvilab.org/), and optionally Python 3.9+ with EasyPQP for database splitting and spectral library generation.
