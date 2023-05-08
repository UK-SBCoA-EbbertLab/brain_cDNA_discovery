#!/bin/bash
#SBATCH --time=23:15:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=align_to_chm13   # Job name
#SBATCH --ntasks=16                 # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --nodes 1
#SBATCH --mem=100G                  # Total memory requested
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err             # Error file for this job.
#SBATCH -o slurm-%j.out             # Output file for this job.
#SBATCH -A cca_mteb223_uksr         # Project allocation account name (REQUIRED)

index=$1
fastq=$2


singularity exec /pscratch/mteb223_uksr/brain_cDNA_discovery/singularity_containers/nanopore.sif \
	minimap2 -t 16 -ax splice \
		-uf \
		$index \
		$fastq > "${fastq}.bam"

samtools sort -@ -12 "${fastq}.bam" -o "${fastq}.sorted.bam"
samtools index "${fastq}.sorted.bam"

samtools flagstat "${fastq}.sorted.bam" -O tsv > "${fastq}.sorted.flagstat.tsv"
samtools idxstats "${fastq}.sorted.bam" > "${fastq}.sorted.idxstat"

#https://groups.google.com/g/rna-star/c/dUKvviBixTQ
samtools view -b -h -F 0x100 -F 0x800 -@ 12 "${fastq}.sorted.bam" > "${fastq}.primary_only.bam"
samtools index "${fastq}.primary_only.bam"
samtools idxstats "${fastq}.primary_only.bam" > "${fastq}.primary_only.idxstat"
samtools flagstat "${fastq}.primary_only.bam" -O tsv > "${fastq}.primary_only.flagstat.tsv"

samtools view -b -q 10 -F 0x100 -F 0x800 -@ 12 "${fastq}.sorted.bam" > "${fastq}.filtered_mapq_primary_alignments.bam"
samtools index "${fastq}.filtered_mapq_primary_alignments.bam"
samtools idxstats "${fastq}.filtered_mapq_primary_alignments.bam" > "${fastq}.filtered_mapq_primary_alignments.idxstat"
samtools flagstat "${fastq}.filtered_mapq_primary_alignments.bam" -O tsv > "${fastq}.filtered_mapq_primary_alignments.flagstat.tsv"


