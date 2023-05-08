# cDNA_pipeline

## Getting Started

### 1) Have a functioning version of Nextflow in your Path.

- Information on how to install NextFlow can be found [here](https://www.nextflow.io/docs/latest/getstarted.html).
          
### 2) Have a functioning version of Singularity on your Path.

- Information on how to install Singularity cna be found [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
          
          
### 3) Clone this github repo using the command below

          git clone https://github.com/UK-SBCoA-EbbertLab/cDNA_pipeline


### 4) Go into the `./workflows/nextflow.config` file and make any necessary changes:

- Alter slurm job manager (or other job manager) parameters to suit your local environment. I don't recommend changing the memory/cpu/time allocated 
for the job manager.
        

          
### 5) Make sure you have all the sequencing files and reference genomes/assemblies files and annotation files you will need to run the pipeline.
          
- ".fastq" -- Nanopore cDNA sequencing files or ".bam" alignment files.

- "sequencing_summary.txt" -- These files are not necessary for execution, but if not available the PycoQC quality control step will be skipped.

- refecence/assembly ".fa" file.

- annotation ".gtf" file is preffered. Only use ".gff3" if using CHM13. Pipeline has an option to handle this, see `Pipeline parameters for STEP 2`.


## Pipeline parameters

          --ont_reads_fastq   <path to fastq sequencing data, can submit multiple at once>
          
          --ont_reads_txt     <path to sequencing summary files, can submit multiple at once. Make sure they follow same naming pattern as fastq files. To ommit enter "None">
          
          --ref               <path to reference/assembly ".fa" file, if using ERCC make sure to concatenate>
  
          --annotation        <path to reference annotation ".gtf" file for GRCh38 or ".gff3" for CHM13. If using GRCh38 and ERCC concatenate the two ".gtf" files>
  
          --out_dir           <name of output directory for pipeline submission. Will appear under "results/<out_dir>" on the directory the pipeline was submitted.
  
          --ercc              <path to ERCC annotation file. Only needed if using CHM13 reference and ERCC. Otherwise set it to "None">
  
          --cdna_kit          <option for pychopper trimming using adapters from the specific cDNA library version, options are "PCS109", "PCS110", "PCS111">
  
          --is_chm13          <set to "True" if using CHM13 and "False" if not. Fixes CHM13 annotation for compatibility with Bambu and converts to ".gtf">
          
          --is_discovery      <logical, if "True" perform transcript discovery and quantification with Bambu, else if "False" perform only
                              quantification based on the given GTF annotation>
                              
          --bambu_track_reads <logical, set to "True" if you want Bambu to keep track of read assignments to transcripts in the output ".RDS" file
                              from Bambu. Set to "False" if you don't need to keep track of read assignments (smaller files). Default: "False">
  


### Submission code used for the analysis contained in this article and repository:

```
nextflow ../main.nf --ont_reads_fq "../../../../../../../scratch/bag222/data/ont_data/09-02-2022_uky_6ad_6ct/*.fastq" \
    --ont_reads_txt "../../../../../../scratch/bag222/data/ont_data/09-02-2022_uky_6ad_6ct/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/" \
    --bambu_track_reads "True" \
    --is_discovery "True" \
    --is_chm13 "False"
 ```
