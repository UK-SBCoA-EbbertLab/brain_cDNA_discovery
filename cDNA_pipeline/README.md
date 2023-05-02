# cDNA_pipeline

## Getting Started

### 1) Have a functioning version of Nextflow in your Path.

- On the MCC it can be found under `/project/mteb223_uksr/sequencing_resources/tools/bin/nextflow`
          
- Make sure to run `module load ccs/java/jdk1.8.0_202` on the MCC load Java. Nextflow needs Java to work. I added this command to my `~/.bash_profile` to make life easier.
          
### 2) Get the singularity image used for this pipeline:

- You can find it on the MCC under: `/project/mteb223_uksr/singularity_files/2022-06-30_cdna_nanopore_pipe.sif`
         
- Alternatively you can build the singularity file from scratch using the ".def" file contained here. I do not recommend this as tools could have been updated and not be compatible with the pipeline anymore.
          
          
### 3) Clone this github repo using the command below

          git clone https://github.com/UK-SBCoA-EbbertLab/cDNA_pipeline


### 4) Put the singularity ".sif" file in the `singularity_container` folder. You can use a softlink if preferred.


### 5) Go into the `workflows/nextflow.config` file and make any necessary changes:

- Alter slurm job manager parameters and singularity file path to suit your needs and PI affiliation. I don't recommend changing the memory/cpu/time allocated 
for the job manager. Particularly, make sure the `--account` slurm job parameters are set to `coa_mteb223_uksr`.
        
- DO NOT change the paths to the files specified by the `params.xxx` variables. Those can be set at the time of executing the pipeline.

          
### 6) Make sure you have all the files, "sequencing_summary.txt" files, reference genomes/assemblies files and annotation files you will need to run the pipeline.
          
- ".fastq" -- cDNA sequencing files.

- "sequencing_summary.txt" --- You can create dummy files if you don't have these available for now. This option will be improved in the future.

- refecence/assembly ".fa" file.

- annotation ".gtf" file is preffered. Only use ".gff3" if using CHM13. Pipeline has an option to handle this.
          


## Pipeline parameters

          --ont_reads_fastq   <path to fastq sequencing data, can submit multiple at once>
          
          --ont_reads_txt     <path to sequencing summary files, can submit multiple at once. Make sure they follow same naming pattern as fastq files. To ommit enter "None">
          
          --ref               <path to reference/assembly ".fa" file, if using ERCC make sure to concatenate>
  
          --annotation        <path to reference annotation ".gtf" file for GRCh38 or ".gff3" for CHM13. If using GRCh38 and ERCC concatenate the two ".gtf" files>
  
          --out_dir           <name of output directory for pipeline submission. Will appear under "results/<out_dir>" on the directory the pipeline was submitted.
  
          --ercc              <path to ERCC annotation file. Only needed if using CHM13 reference and ERCC. Otherwise set it to "None">
  
          --cdna_kit          <option for pychopper trimming using adapters from the specific cDNA library version, options are "PCS109", "PCS110", "PCS111">
  
          --is_chm13          <set to "True" if using CHM13 and "False" if not. Fixes CHM13 annotation for compatibility with Bambu and converts to ".gtf">
  
  


## Examples of the submission of the pipeline can be seen found here `workflow/test_workflows/`, see them below:


### CHM13 without ERCCs

          nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.fastq" \
              --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.txt" \
              --ref "../../references/chm13v2.0.fa" \
              --annotation "../../references/CHM13.v2.0.gff3" \
              --ercc "None" \
              --out_dir "./CHM13_test/" \
              --cdna_kit "PCS111" \
              --is_chm13 "True"  -resume

### CHM13 with ERCCs

          nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.fastq" \
              --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.txt" \
              --ref "../../references/chm13v2.0_ERCC.fa" \
              --annotation "../../references/CHM13.v2.0.gff3" \
              --ercc "../../references/ERCC92.gtf" \
              --out_dir "./CHM13_ERCC_test/" \
              --cdna_kit "PCS111" \
              --is_chm13 "True"  -resume
    
#### Notice that for CHM13 you need to concatenate CHM13 and the ERCC reference prior to submitting the pipeline, but the annotations are entered separately and concatenated by the program itself after converting CHM13 annotation to ".gtf" format.

### GRCh38 without ERCCs

          nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.fastq" \
              --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.txt" \
              --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
              --annotation "../../references/Homo_sapiens.GRCh38.106.gtf" \
              --out_dir "./GRCh38_test/" \
              --ercc "None" \
              --cdna_kit "PCS111" \
              --is_chm13 "False"  -resume


### GRCh38 with ERCCs

          nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.fastq" \
              --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.txt" \
              --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
              --annotation "../../references/Homo_sapiens.GRCh38.106_ERCC.gtf" \
              --out_dir "./GRCh38_ERCC_test/" \
              --ercc "None" \
              --cdna_kit "PCS111" \
              --is_chm13 "False"  -resume
    
#### Notice that for GRCh38 the `--ercc` is always set to "None" as there the user can easily concatenate both the GRCh38 reference and the annotation to the ERCC reference and annotation.          
