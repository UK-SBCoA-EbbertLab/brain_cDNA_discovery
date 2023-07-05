# Article Title

## Link to pre-print:

Link

### This repository contains all code and documentation used for the analysis contained in the article above.

## Repository structure:

`article_analysis` - Scripts used for data analysis and figure generation for article publication. Uses data from output from `cDNA_pipeline`, and `proteomics_pipeline`, Some figures were created using the scripts in `website`


`cDNA_pipeline` - In house NextFlow pipeline optimized for analysis of Oxford Nanopore PCR Amplified cDNA sequencing data.


`proteomics_pipeline` - Fragpipe pipeline used to validate new transcripts at the protein level using public Mass Spec data.

`singularity_containers` - Directory with container definition files and pull commands. With the exception of the Fragpipe pipeline (`proteomics_pipeline`) and the Rshiny web app (`website`), all the software used in this GitHub repository is in these singularity containers.

`website` - Contains Rshiny app scripts that allows users to perform gene queries and visualize RNA isoform expression from the data used in this publication.
URL: 

## Data availability

Raw nanopore sequencing fastq files generated in this study are available at: XXX

Final output files and annotations/references used in this study are available at: XXX

Public proteomics (Mass spec) data from cell-lines used in this experiment are available [here](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024364). For more information about this data see: https://pubmed.ncbi.nlm.nih.gov/36959352/

Public proteomics (Mass spec) data from brain are available [here](https://www.synapse.org/#!Synapse:syn25006611/wiki/608683). For more information about this data see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8825285/

## More information

Each directory within this GitHub repository contains documentation for the analysis performed in that directory.
If you have any questions please submit and issue.
