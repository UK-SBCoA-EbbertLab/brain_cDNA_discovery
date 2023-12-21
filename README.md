# Using deep long-read RNAseq in Alzheimer’s disease brain to assess medical relevance of RNA isoform diversity

## Link to pre-print: https://www.biorxiv.org/content/10.1101/2023.08.06.552162v2

### This repository contains all code and documentation used for the analysis contained in the article above.

## Repository structure:

`article_analysis` - Scripts used for data analysis and figure generation for article publication. Uses data from output from `illumina_pipeline` and `nanopore_pipeline`. Some figures were created using the scripts in `website`

`illumina_pipeline` - In house NextFlow pipeline optimezed for analysis of Illumina paired-end short-read sequencing data.

`nanopore_pipeline` - In house NextFlow pipeline optimized for analysis of Oxford Nanopore PCR Amplified cDNA sequencing data.


`proteomics` - Analysis pipeline to validate new transcripts at the protein level using publicly available Mass Spec data. Also explains downstream analysis steps and
contains custom script used for downstream analysis and figure generation.

`singularity_containers` - Directory with container definition files and pull commands. With the exception of the Fragpipe pipeline (`proteomics_pipeline`) and the Rshiny web app (`website`), all the software used in this GitHub repository is in these singularity containers.

`website` - Contains Rshiny app scripts that allows users to perform gene queries and visualize RNA isoform expression from the data used in this publication.
Access website [here](https://ebbertlab.com/brain_rna_isoform_seq.html)

## Data availability

Raw nanopore sequencing fastq files generated in this study are available [here](https://www.synapse.org/#!Synapse:syn52047893/wiki/622953). Also available through NIH SRA (Accession number: SRP456327)

Proteomics (Mass spec) data from cell-lines used in this experiment are publicly available [here](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024364). For more information about this data see: https://pubmed.ncbi.nlm.nih.gov/36959352/

Proteomics (Mass spec) data from round 2 of the ROSMAP TMT brain Proteomics are puclicly available [here](https://www.synapse.org/#!Synapse:syn17015098). For more information about this data see: https://www.nature.com/articles/s41597-020-00650-8

Final output files from transcriptomics/RNAseq and proteomics analysis and annotations/references used in this study are available [here](https://doi.org/10.5281/zenodo.8180677)

GTEx long-read RNAseq data used for validation of our study results is available [here](https://anvil.terra.bio/#workspaces/anvil-datastorage/AnVIL_GTEx_V9_hg38)

ROSMAP short-read RNAseq data used for validation of our study results is available [here](https://www.synapse.org/#!Synapse:syn21589959)


## More information

Each directory within this GitHub repository contains documentation for the analysis performed in that directory.
If you have any questions please submit and issue.
