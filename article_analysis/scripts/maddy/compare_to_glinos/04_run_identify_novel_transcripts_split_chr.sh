#!/bin/bash
#SBATCH --time=15:15:00                                         # Time limit for the job (REQUIRED).
#SBATCH --job-name=GTF_comparisons      # Job name
#SBATCH --ntasks=1                                              # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=100G                                                # Total memory requested
#SBATCH --partition=normal                                      # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err                                         # Error file for this job.
#SBATCH -o slurm-%j.out                                         # Output file for this job.
#SBATCH -A coa_mteb223_uksr                                     # Project allocation account name (REQUIRED)


### OURS vs GLINOS
#our_data="compare_to_glinos/ad_vs_ct_pilot_annotations_just_novel.gtf"
#their_data="compare_to_glinos/GTEx_annotations_just_novel.gtf"
#output_file_prefix="Ours_vs_GTEx_condensed_threshold_6_9_23"
#novel_column_alt_value="GTEX"

### OURS vs ENSEMBL
#our_data="extended_annotation_grch38_only_novel.sorted.gtf"
#their_data="Homo_sapiens.GRCh38.105.simplified.gtf"
#output_file_prefix="Ours_vs_Ensembl_condensed_threshold"
#novel_column_alt_value="ANNOTATED"

### GLINOS vs OURS
our_data="compare_to_glinos/GTEx_annotations_just_novel.gtf"
their_data="compare_to_glinos/ad_vs_ct_pilot_annotations_just_novel.gtf"
output_file_prefix="Glinos_vs_Ours_condensed_threshold_6_9_23"
novel_column_alt_value="OURS"


### GLINOS vs ENSEMBL
#our_data="flair_filter_transcripts_no_chr.sorted.gtf"
#their_data="Homo_sapiens.GRCh38.105.simplified.gtf"
#output_file_prefix="Glinos_vs_Ensembl_condensed_threshold"
#novel_column_alt_value="ANNOTATED"


for seqname in {1..22}
do
	awk -v num=$seqname -v FS='\t' -v OFS='\t' '$1 == num' ${our_data} > "compare_to_glinos/compare_novel_transcript_annotations/.${output_file_prefix}.our_data.chr_${seqname}.tmp"
	awk -v num=$seqname -v FS='\t' -v OFS='\t' '$1 == num' ${their_data} > "compare_to_glinos/compare_novel_transcript_annotations/.${output_file_prefix}.their_data.chr_${seqname}.tmp"
	sbatch -n 1 -t 04:00:00 -A coa_mteb223_uksr -p normal -J comparision_${seqname} --mem=100G --wrap=" \
	time singularity exec /project/mteb223_uksr/singularity_files/transcriptome_long_reads_GTEx_2022_12_14.sif \
	Rscript 04_identify_novel_transcripts.R \
		compare_to_glinos/compare_novel_transcript_annotations/.${output_file_prefix}.our_data.chr_${seqname}.tmp \
		compare_to_glinos/compare_novel_transcript_annotations/.${output_file_prefix}.their_data.chr_${seqname}.tmp \
		compare_to_glinos/compare_novel_transcript_annotations/${output_file_prefix}_chr_${seqname}.tsv \
		$novel_column_alt_value"
done


awk -v FS='\t' -v OFS='\t' !'/^[0-9]/' ${our_data} > "compare_to_glinos/compare_novel_transcript_annotations/.${output_file_prefix}.our_data.chr_non_numeric.tmp"
awk -v FS='\t' -v OFS='\t' !'/^[0-9]/' ${their_data} > "compare_to_glinos/compare_novel_transcript_annotations/.${output_file_prefix}.their_data.chr_non_numeric.tmp"
sbatch -n 1 -t 04:00:00 -A coa_mteb223_uksr -p normal -J comparision_non_numeric --mem=100G --wrap=" \
time singularity exec /project/mteb223_uksr/singularity_files/transcriptome_long_reads_GTEx_2022_12_14.sif \
	Rscript 04_identify_novel_transcripts.R \
	compare_to_glinos/compare_novel_transcript_annotations/.${output_file_prefix}.our_data.chr_non_numeric.tmp \
	compare_to_glinos/compare_novel_transcript_annotations/.${output_file_prefix}.their_data.chr_non_numeric.tmp \
	compare_to_glinos/compare_novel_transcript_annotations/${output_file_prefix}_chr_non_numeric.tsv \
	$novel_column_alt_value"


