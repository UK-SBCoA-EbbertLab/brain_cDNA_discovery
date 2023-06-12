#!/bin/bash

### Take the output from 04_identify_novel_transcripts.R and get the number of transcripts at each cutoff
### Then create a line graph with all the numbers
echo -e "File_name\tDiff_1bp\tDiff_5bp\tDiff_10bp\tDiff_20bp\tDiff_50bp\tDiff_100bp\tDiff_200bp" > ../summary_tibble_file.tsv
echo -e "File_name\tDiff_1bp\tDiff_5bp\tDiff_10bp\tDiff_20bp\tDiff_50bp\tDiff_100bp\tDiff_200bp" > ../summary_tibble_file_normalized.tsv
#TODO fix these paths
for file in ../glinos_vs_ensembl/Glinos_vs_Ensembl_condensed_threshold.tsv ../glinos_vs_ours/Glinos_vs_Ours_condensed_threshold.all.tsv ../ours_vs_ensembl/Ours_vs_Ensembl_condensed_threshold.tsv ../ours_vs_glinos/Ours_vs_Glinos_condensed_threshold.all.tsv
do
	echo $file
	bp_1=$(awk -v FS='\t' -v OFS='\t' '$8 != "NOVEL" {print}' $file | wc -l)
	bp_5=$(awk -v FS='\t' -v OFS='\t' '$9 != "NOVEL" {print}' $file | wc -l)
	bp_10=$(awk -v FS='\t' -v OFS='\t' '$10 != "NOVEL" {print}' $file | wc -l)
	bp_20=$(awk -v FS='\t' -v OFS='\t' '$11 != "NOVEL" {print}' $file | wc -l)
	bp_50=$(awk -v FS='\t' -v OFS='\t' '$12 != "NOVEL" {print}' $file | wc -l)
	bp_100=$(awk -v FS='\t' -v OFS='\t' '$13 != "NOVEL" {print}' $file | wc -l)
	bp_200=$(awk -v FS='\t' -v OFS='\t' '$14 != "NOVEL" {print}' $file | wc -l)
	file_length=$(wc -l <$file)

	bp_1_n=$(awk -v var1=$bp_1 -v var2=$file_length 'BEGIN { print  ( var1 / var2 * 100 ) }')
	bp_5_n=$(awk -v var1=$bp_5 -v var2=$file_length 'BEGIN { print  ( var1 / var2 * 100 ) }')
	bp_10_n=$(awk -v var1=$bp_10 -v var2=$file_length 'BEGIN { print  ( var1 / var2 * 100 ) }')
	bp_20_n=$(awk -v var1=$bp_20 -v var2=$file_length 'BEGIN { print  ( var1 / var2 * 100 ) }')
	bp_50_n=$(awk -v var1=$bp_50 -v var2=$file_length 'BEGIN { print  ( var1 / var2 * 100 ) }')
	bp_100_n=$(awk -v var1=$bp_100 -v var2=$file_length 'BEGIN { print  ( var1 / var2 * 100 ) }')
	bp_200_n=$(awk -v var1=$bp_200 -v var2=$file_length 'BEGIN { print  ( var1 / var2 * 100 ) }')

	echo -e "${file}\t${bp_1}\t${bp_5}\t${bp_10}\t${bp_20}\t${bp_50}\t${bp_100}\t${bp_200}" >> ../summary_tibble_file.tsv
	echo -e "${file}\t${bp_1_n}\t${bp_5_n}\t${bp_10_n}\t${bp_20_n}\t${bp_50_n}\t${bp_100_n}\t${bp_200_n}" >> ../summary_tibble_file_normalized.tsv
done

singularity exec /project/mteb223_uksr/singularity_files/transcriptome_long_reads_GTEx.sif Rscript 05a_graph_summarize_novel.R ../summary_tibble_file.tsv ../summary_tibble_file_normalized.tsv
