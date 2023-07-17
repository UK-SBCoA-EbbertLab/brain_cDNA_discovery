# proteomics\_pipeline

## Getting Started

### 1) Install Fragpipe and necessary dependencies[^1].

 - Information on how to install Fragpipe can be found [here](https://fragpipe.nesvilab.org/).

### 2) Download public proteomics (mass spec) data.

 - Data from cell-lines can be found [here](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024364).
	
 - Data from brain can be found [here]([200~https://www.synapse.org/#!Synapse:syn17015098).

### 3) Clone this repository using the command below:

```
	git clone https://github.com/UK-SBCoA-EbbertLab/brain_cDNA_discovery
```

### 4) Modify `fragpipe.workflow`, and `fragpipe_submission_script`:

 - Use X11 Forwarding to load the workflow into the Fragpipe GUI and make necessary changes.
 
 - Modify the submission script to include the correct paths and working directory. 

 - If using a workload manager then give the job enough time to run (36 hours should be more than enough).

## Downstream Analysis

### Cell-lines

	1) Open the `combined_peptide.tsv` in Excel.
	
	2) Filter the "Protein" column for IDs containing "Bambu" and filter the "Mapped Genes" and "Mapped Proteins" columns for blank entries.

	3) Run `count_peptides.sh` with the remaining peptide sequences and `peptide_database_uniprot_header_run_MEDIAN_cpm_greater_than_one.fa` to get the number of times each peptide sequence shows up in the database.

	4) Filter out any peptide sequences that had a count > 1.
	
	5) Create a column to sum each spectral count entry for each peptide from each sample.

	6) Group by "Protein".

	7) Filter the total spectral count to only show values > 5.

	8) Use the remaining protein IDs to get corresponding list of peptide sequences.

	9) Run `create_fasta.sh` with the list of protein IDs and peptide sequences to produce a fasta file.

### Brain
	
	Since there are 14 batches for this data, downstream analysis is slightly different.

	1) Modify `merge_peptides.py` with the correct directory for the peptide.tsv files.

	2) Open `merge_peptides.py` as a jupyter notebook, run,  and export the column of peptide sequences as a csv file.

	3) Run `count_peptides` with the peptide sequences csv file that was just created and `peptide_database_uniprot_header_run_MEDIAN_cpm_greater_than_one.fa` to get the number or times each peptide sequence shows up in the database.

	4) Filter out any peptide sequence that had a count > 1.

	5) Filter "Spectral Count" for values > 5.
	
	6) Run `create_fasta.sh` with  the remaining protein IDs and corresponding peptide sequences to produce a fasta file.

[^1]:This will include Java 9+, [MSFragger](https://msfragger.nesvilab.org/), [Philosopher](https://philosopher.nesvilab.org/), [IonQuant](http://ionquant.nesvilab.org/), and optionally Python 3.9+ with EasyPQP for database splitting and spectral library generation.
