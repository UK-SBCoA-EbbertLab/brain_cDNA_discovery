# nanopore_pipeline

## Getting Started

### 1) Have a functioning version of Nextflow in your Path.

- Information on how to install NextFlow can be found [here](https://www.nextflow.io/docs/latest/getstarted.html).
          
### 2) Have a functioning version of Singularity on your Path.

- Information on how to install Singularity cna be found [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
          
          
### 3) Clone this github repo using the command below

          git clone https://github.com/UK-SBCoA-EbbertLab/brain_cDNA_discovery


### 4) Go into the `./workflows/nextflow.config` file and make any necessary changes:

- Alter slurm job manager (or other job manager) parameters to suit your local environment. I don't recommend changing the memory/cpu/time allocated 
for the job manager. Also, do not change the singularity parameters. Those commands will automatically pull the necessary containers to run the analysis.
        

          
### 5) Make sure you have all the sequencing files and reference genomes/assemblies files and annotation files you will need to run the pipeline.
          
- "R1.fastq.gz" and "R2.fastq.gz" -- Paired-end illumina RNAseq fastq files (gzipped is preferred).

- refecence/assembly ".fa" file.

- annotation ".gtf" file is preffered.

- housekeeping gene/transcripts "hg38.HouseKeepingGenes.bed"

- transcriptome ".fa" file. Can be made using [gffread](https://github.com/gpertea/gffread) to process the assembly ".fa" file with the annotation ".gtf" file using this command:
`gffread -w transcripts.fa -g reference_genome.fa annotation.gtf` 

See home directory for this repository to find data needed to reproduce the analysis done in this article.


## Pipeline parameters

          --illumina_data   <path to illumina paired-end RNAseq ".fastq" data. 
                             To make sure paired-end files from same sample are processed together path should be specified ub this format: /path/"*_R{1,2}*.fastq.gz">
                             
          --ref               <path to reference/assembly ".fa" file, if using ERCC make sure to concatenate>
    
          --out_dir           <name of output directory for pipeline submission. Will appear under "results/<out_dir>" on the directory the pipeline was submitted.
    
          --overhang          <This should be the size of your paired-end reads - 1. So for 150bp paired end reads, this parameter should be set to 149>
            
          --transcriptome      <path to reference transcriptome ".fa" file>
                              
          --bam                <path to aligned ".bam" files. If you have files that are already aligned you can use this parameter to filter reads and only keep reads that uniquely align one transcript
                              in your transcriptome file (MAPQ>255). If you specify this parameter only need to specify a the --transcriptome parameter as a complement, all other parameters will be ignored>
  


### Submission code used for initial ROSMAP analysis

```
nextflow ../main.nf --illumina_data "/../../../../../scratch/bag222/data/ROSMAP_illumina_RNAseq/DorsoLateralFrontalCortex/*_R{1,2}*.fastq.gz" \
    --ref "../../../nanopore_pipeline/references/Homo_sapiens.GRCh38_ERCC.fa" \
    --transcriptome "../../../article_analysis/data/raw/nextflow_pipeline_output/transcriptome/transcriptome.fa" \
    --annotation "../../../article_analysis/data/raw/nextflow_pipeline_output/bambu_discovery/extended_annotations.gtf" \
    --housekeeping "../../../cDNA_pipeline/references/hg38.HouseKeepingGenes.bed" \
    --out_dir "./ROSMAP_illumina_DorsoLateralPreFrontalCortex_with_our_extended_annotation/" \
    --overhang "149" -resume
 ```

### Submission code used to only keep reads that uniquely align one transcript in ROSMAP data

```
nextflow ../main.nf --bam "./results/ROSMAP_illumina_DorsoLateralPreFrontalCortex_with_our_extended_annotation/STAR/*toTranscriptome.out.bam" \
    --transcriptome "../../../article_analysis/data/raw/nextflow_pipeline_output/transcriptome/transcriptome.fa" \
    --out_dir "./ROSMAP_illumina_DorsoLateralPreFrontalCortex_UNIQUE/" -resume
```

## Directory structure:

`modules` - NextFlow module ".nf" files for different steps in the pipeline.

`subworkflows` - NextFlow subworkflow ".nf" files for different executions of the pipeline.

`workflow` - Contains main NextFlow workflow, configuration files, submission scripts, and custom R/Python scripts used in the pipeline.

`CITATION.md` - Citations to tools used in the NextFlow pipeline.

`LICENSE` - License for NextFlow pipeline.


## Pipeline overview (submission used in the article)


  1) Adapter trimming and quality control of fastq files with `trim_galore`.
  2) Alignment to reference transcriptome using `STAR`.
  3) Transcript quantification with `salmon`.

## Alternative submission
  1) Filtering of aligned reads that are not unique (only keeps reads with MAPQ = 255) with `samtools`
  2) Transcript quantification with `salmon`

## More information:

We obtained publicly available ROSMAP Illumina 150bp paired-end RNAseq data from 251 brain samples (Brodmann area 9/46). Samples were processed using this `illumina_pipeline` NextFlow pipeline. More detailed information about the samples, sequencing protocol, and analysis pipeline can be found in the article publication associated with this project.

