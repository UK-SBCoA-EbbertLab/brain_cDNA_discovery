#!/bin/bash

#sort -k1,1 -k4,4n Homo_sapiens.GRCh38.105.simplified.gtf > Homo_sapiens.GRCh38.105.simplified.sorted.gtf
#sort -k1,1 -k4,4n extended_annotation_grch38_only_novel.gtf > extended_annotation_grch38_only_novel.sorted.gtf
#sort -k1,1 -k4,4n flair_filter_transcripts_no_chr.gtf > flair_filter_transcripts_no_chr.sorted.gtf

#awk -v FS='\t' -v OFS='\t' '$3 == "exon" {print}' Homo_sapiens.GRCh38.105.simplified.sorted.gtf > Homo_sapiens.GRCh38.105.simplified.sorted.just_exons.gtf
#awk -v FS='\t' -v OFS='\t' '$3 == "exon" {print}' extended_annotation_grch38_only_novel.sorted.gtf > extended_annotation_grch38_only_novel.sorted.just_exons.gtf
#awk -v FS='\t' -v OFS='\t' '$3 == "exon" {print}' flair_filter_transcripts_no_chr.sorted.gtf > flair_filter_transcripts_no_chr.sorted.just_exons.gtf
#awk -v FS='\t' -v OFS='\t' '$3 == "transcript" {print}' Homo_sapiens.GRCh38.105.simplified.sorted.gtf > Homo_sapiens.GRCh38.105.simplified.sorted.just_transcripts.gtf
#awk -v FS='\t' -v OFS='\t' '$3 == "transcript" {print}' extended_annotation_grch38_only_novel.sorted.gtf > extended_annotation_grch38_only_novel.sorted.just_transcripts.gtf
#awk -v FS='\t' -v OFS='\t' '$3 == "transcript" {print}' flair_filter_transcripts_no_chr.sorted.gtf > flair_filter_transcripts_no_chr.sorted.just_transcripts.gtf


outfile_prefix=("ours_vs_ensembl" "ours_vs_glinos" "glinos_vs_ensembl" "glinos_vs_ours")
pattern_regex=('/gene.[0-9]/' '/gene.[0-9]/' '/\"chr[0-9]+:[0-9]+\"/' '/\"chr[0-9]+:[0-9]+\"/')
a_gtf_transcript=("../gtfs/extended_annotation_grch38_only_novel.sorted.just_transcripts.gtf" "../gtfs/extended_annotation_grch38_only_novel.sorted.just_transcripts.gtf" "../gtfs/flair_filter_transcripts_no_chr.sorted.just_transcripts.gtf" "../gtfs/flair_filter_transcripts_no_chr.sorted.just_transcripts.gtf")
b_gtf_transcript=("../gtfs/Homo_sapiens.GRCh38.105.simplified.sorted.just_transcripts.gtf" "../gtfs/flair_filter_transcripts_no_chr.sorted.just_transcripts.gtf" "../gtfs/Homo_sapiens.GRCh38.105.simplified.sorted.just_transcripts.gtf" "../gtfs/extended_annotation_grch38_only_novel.sorted.just_transcripts.gtf")
a_gtf_exon=("../gtfs/extended_annotation_grch38_only_novel.sorted.just_exons.gtf" "../gtfs/extended_annotation_grch38_only_novel.sorted.just_exons.gtf" "../gtfs/flair_filter_transcripts_no_chr.sorted.just_exons.gtf" "../gtfs/flair_filter_transcripts_no_chr.sorted.just_exons.gtf")
b_gtf_exon=("../gtfs/Homo_sapiens.GRCh38.105.simplified.sorted.just_exons.gtf" "../gtfs/flair_filter_transcripts_no_chr.sorted.just_exons.gtf" "../gtfs/Homo_sapiens.GRCh38.105.simplified.sorted.just_exons.gtf" "../gtfs/extended_annotation_grch38_only_novel.sorted.just_exons.gtf")

for i in {0..3};
do
	for j in None 0.5 0.9;
	do
		bash compare_novel_to_annotated.sh ${a_gtf_transcript[i]} ${b_gtf_transcript[i]} ${outfile_prefix[i]} "just_transcripts" ${pattern_regex[i]} $j
		bash compare_novel_to_annotated.sh ${a_gtf_exon[i]} ${b_gtf_exon[i]} ${outfile_prefix[i]} "just_exons" ${pattern_regex[i]} $j
	done

done

