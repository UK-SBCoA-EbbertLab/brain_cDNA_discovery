# Article Title

## Link to pre-print:

Link

### This repository contains all code and documentation used for the analysis contained in the article above.

## Repository structure:

`article_analysis` - Scripts used for data analysis and figure generation for article publication. Uses data from output from `cDNA_pipeline`, and `proteomics_pipeline`, Some figures were created using the scripts in `website`


`cDNA_pipeline` - In house NextFlow pipeline optimized for analysis of Oxford Nanopore PCR Amplified cDNA sequencing data.


`proteomics` - Analysis pipeline to validate new transcripts at the protein level using publicly available Mass Spec data.

`singularity_containers` - Directory with container definition files and pull commands. With the exception of the Fragpipe pipeline (`proteomics_pipeline`) and the Rshiny web app (`website`), all the software used in this GitHub repository is in these singularity containers.

`website` - Contains Rshiny app scripts that allows users to perform gene queries and visualize RNA isoform expression from the data used in this publication.
Access website [here](https://ebbertlab.shinyapps.io/Transcripts_and_counts/)

## Data availability

Raw nanopore sequencing fastq files generated in this study are available [here](https://www.synapse.org/ebbert_lab_brain_long_read_cDNA_discovery_project)

Proteomics (Mass spec) data from cell-lines used in this experiment are publicly available [here](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024364). For more information about this data see: https://pubmed.ncbi.nlm.nih.gov/36959352/

Proteomics (Mass spec) data from round 2 of the ROSMAP TMT brain Proteomics are puclicly available [here](https://www.synapse.org/#!Synapse:syn17015098). For more information about this data see: https://www.nature.com/articles/s41597-020-00650-8

Final output files from transcriptomics/RNAseq and proteomics analysis and annotations/references used in this study are available at: <zenodo_link>


## More information

Each directory within this GitHub repository contains documentation for the analysis performed in that directory.
If you have any questions please submit and issue.
