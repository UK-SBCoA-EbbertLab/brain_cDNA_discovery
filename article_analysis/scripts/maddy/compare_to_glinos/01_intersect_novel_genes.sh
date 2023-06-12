#!/bin/bash

### The purpose of this script is to assess the overlap of our novel transcripts with new geneids with glinos transcripts with new geneids.
### This is looking purely at the whole transcript level and not at individual exons. 
### A 90% overlap of either our transcript or glinos' transcript is required.

our_gtf='compare_to_glinos/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/bambu_discovery/extended_annotations.gtf'
gtex_gtf='compare_to_glinos/GTEx_GRCh38_version_107_mapq_10_track_reads/bambu_discovery/extended_annotations.gtf'
our_gtf_sorted='compare_to_glinos/ad_vs_ct_pilot_study_feb_2023_GRCh38-107.extended_annotation.sorted.gtf'
gtex_gtf_sorted='compare_to_glinos/GTEx_GRCh38-107.extended_annotation.sorted.gtf'

our_gtf_sorted_novel_genes='compare_to_glinos/ad_vs_ct_pilot_study_feb_2023_GRCh38-107.extended_annotation.sorted.only_novel_genesgtf'
gtex_gtf_sorted_novel_genes='compare_to_glinos/GTEx_GRCh38-107.extended_annotation.only_novel_genes.sorted.gtf'

our_stranded_95_overlap='compare_to_glinos/our_95_gtf_overlap_stranded.gtf'
our_unstranded_95_overlap='compare_to_glinos/our_95_gtf_overlap_not_stranded.gtf'

gtex_stranded_95_overlap='compare_to_glinos/gtex_95_gtf_overlap_stranded.gtf'
gtex_unstranded_95_overlap='compare_to_glinos/gtex_95_gtf_overlap_not_stranded.gtf'


######### PREP THE GTFS ################################################

# Ensure the gtf's are sorted
sort -k1,1 -k4,4n $our_gtf > $our_gtf_sorted
sort -k1,1 -k4,4n $gtex_gtf > $gtex_gtf_sorted

## Grab only the entries whose feature states 'transcript' and grab the entries without an ENSMBL id
awk -F "\t" '$3 == "transcript"' ${our_gtf_sorted} | awk -F "; " '$1 !~ "ENSG"' > ${our_gtf_sorted_novel_genes}
awk -F "\t" '$3 == "transcript"' ${gtex_gtf_sorted} | awk -F "; " '$1 !~ "ENSG"' > ${gtex_gtf_sorted_novel_genes}

###########################################################################



######## ASSESS WRT STRANDEDNESS, REPORT OUR TRANSCRIPTS #################
 
# interesct the gtfs that have only novel genes. Keep full features from our data that overlap
# Overlap must be at least 90% for one of the features
# THIS IS STRANDED
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
       -wa \
       -F 0.95 \
       -f 0.95 \
       -e \
       -u \
       -s \
       -a ${our_gtf_sorted_novel_genes} \
       -b ${gtex_gtf_sorted_novel_genes} \
       > ${our_stranded_95_overlap}

# subtract out the overlap and get just our uniquely novel transcripts
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools subtract \
	-s \
	-a ../gtfs/extended_annotation_grch38_only_novel_genes.sorted.gtf \
	-b ../overlap_assessment/just_novel_genes_our_overlap_stranded.gtf \
	> ../overlap_assessment/our_novel_genes_stranded.gtf

# upon merging, the file becomes a bed file format
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools merge \
	-s \
	-i ../overlap_assessment/our_novel_genes_stranded.gtf \
	> ../overlap_assessment/our_novel_genes_stranded.merged.bed

# This creates an additional column of gap lengths to be graphed later on
python3 01a_create_gap_lengths_column.py ../overlap_assessment/our_novel_genes_stranded.merged.bed

######################################################################




######### ASSESS WITHOUT STRANDEDNESS, REPORT OUR TRANSCRIPTS ################

# interesct the gtfs that have only novel genes. Keep full features from our data that overlap
# Overlap must be at least 90% for one of the features
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
       -wa \
       -F 0.9 \
       -f 0.9 \
       -e \
       -u \
       -a ../gtfs/extended_annotation_grch38_only_novel_genes.sorted.gtf \
       -b ../gtfs/flair_filter_transcripts_no_chr_novel_genes.sorted.gtf \
       > ../overlap_assessment/just_novel_genes_our_overlap.gtf

# subtract out the overlap and get just our uniquely novel transcripts
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools subtract \
        -a ../gtfs/extended_annotation_grch38_only_novel_genes.sorted.gtf \
        -b ../overlap_assessment/just_novel_genes_our_overlap.gtf \
        > ../overlap_assessment/our_novel_genes.gtf


# upon merging, the file becomes a bed file format
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools merge \
        -i ../overlap_assessment/our_novel_genes.gtf \
        > ../overlap_assessment/our_novel_genes.merged.bed

# This creates an additional column of gap lengths to be graphed later on
python3 01a_create_gap_lengths_column.py ../overlap_assessment/our_novel_genes.merged.bed

#############################################################################




######## ASSESS WRT STRANDEDNESS, REPORT OUR TRANSCRIPTS AGAINST ENSEMBL #################
 
# interesct the gtfs that have only novel genes. Keep full features from our data that overlap
# Overlap must be at least 90% for one of the features
# THIS IS STRANDED
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
       -wa \
       -F 0.9 \
       -f 0.9 \
       -e \
       -u \
       -s \
       -a ../gtfs/extended_annotation_grch38_only_novel_genes.sorted.gtf \
       -b ../gtfs/Homo_sapiens.GRCh38.105.simplified.sorted.gtf\
       > ../overlap_assessment/just_novel_genes_our_overlap_stranded_against_ensembl.gtf

# subtract out the overlap and get just our uniquely novel transcripts
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools subtract \
	-s \
	-a ../gtfs/extended_annotation_grch38_only_novel_genes.sorted.gtf \
	-b ../overlap_assessment/just_novel_genes_our_overlap_stranded_against_ensembl.gtf \
	> ../overlap_assessment/our_novel_genes_stranded_against_ensembl.gtf

# upon merging, the file becomes a bed file format
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools merge \
	-s \
	-i ../overlap_assessment/our_novel_genes_stranded_against_ensembl.gtf \
	> ../overlap_assessment/our_novel_genes_stranded_against_ensembl.merged.bed

######################################################################




######### ASSESS WITHOUT STRANDEDNESS, REPORT OUR TRANSCRIPTS AGAINST ENSEMBL ################

# interesct the gtfs that have only novel genes. Keep full features from our data that overlap
# Overlap must be at least 90% for one of the features
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
       -wa \
       -F 0.9 \
       -f 0.9 \
       -e \
       -u \
       -a ../gtfs/extended_annotation_grch38_only_novel_genes.sorted.gtf \
       -b ../gtfs/Homo_sapiens.GRCh38.105.simplified.sorted.gtf \
       > ../overlap_assessment/just_novel_genes_our_overlap_against_ensembl.gtf

# subtract out the overlap and get just our uniquely novel transcripts
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools subtract \
        -a ../gtfs/extended_annotation_grch38_only_novel_genes.sorted.gtf \
        -b ../overlap_assessment/just_novel_genes_our_overlap_against_ensembl.gtf \
        > ../overlap_assessment/our_novel_genes_against_ensembl.gtf


# upon merging, the file becomes a bed file format
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools merge \
        -i ../overlap_assessment/our_novel_genes_against_ensembl.gtf \
        > ../overlap_assessment/our_novel_genes_against_ensembl.merged.bed

##########################################################################




######### ASSESS GLINOS WRT STRANDEDNESS, AGAINST ENSEMBL ##############

# interesct the gtfs that have only novel genes. Keep full features from glinos data that overlap
# Overlap must be at least 90% for one of the features
# THIS IS STRANDED
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
       -wa \
       -F 0.9 \
       -f 0.9 \
       -s \
       -e \
       -u \
       -a ../gtfs/flair_filter_transcripts_no_chr_novel_genes.sorted.gtf \
       -b ../gtfs/Homo_sapiens.GRCh38.105.simplified.sorted.gtf \
       > ../overlap_assessment/just_novel_genes_their_overlap_stranded_against_ensembl.gtf

#subtract out the overlap and get just their uniquely novel transcripts
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools subtract \
	-s \
        -a ../gtfs/flair_filter_transcripts_no_chr_novel_genes.sorted.gtf \
        -b ../overlap_assessment/just_novel_genes_their_overlap_stranded_against_ensembl.gtf \
        > ../overlap_assessment/their_novel_genes_stranded_against_ensembl.gtf


# upon merging, the file becomes a bed file format
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools merge \
	-s \
        -i ../overlap_assessment/their_novel_genes_stranded_against_ensembl.gtf \
        > ../overlap_assessment/their_novel_genes_stranded_against_ensembl.merged.bed

##########################################################################




########## ASSESS GLINOS WITHOUT STRANDEDNESS, AGAINST ENSEMBL ###########
 
# interesct the gtfs that have only novel genes. Keep full features from glinos data that overlap
# Overlap must be at least 90% for one of the features
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
       -wa \
       -F 0.9 \
       -f 0.9 \
       -e \
       -u \
       -a ../gtfs/flair_filter_transcripts_no_chr_novel_genes.sorted.gtf \
       -b ../gtfs/Homo_sapiens.GRCh38.105.simplified.sorted.gtf \
       > ../overlap_assessment/just_novel_genes_their_overlap_against_ensembl.gtf

#subtract out the overlap and get just their uniquely novel transcripts
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools subtract \
	-a ../gtfs/flair_filter_transcripts_no_chr_novel_genes.sorted.gtf \
	-b ../overlap_assessment/just_novel_genes_their_overlap_against_ensembl.gtf \
	> ../overlap_assessment/their_novel_genes_against_ensembl.gtf


# upon merging, the file becomes a bed file format
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools merge \
	-i ../overlap_assessment/their_novel_genes_against_ensembl.gtf \
	> ../overlap_assessment/their_novel_genes_against_ensembl.merged.bed

#########################################################################




######### ASSESS WRT STRANDEDNESS, REPORT GLINOS TRANSCRIPTS ##############

# interesct the gtfs that have only novel genes. Keep full features from glinos data that overlap
# Overlap must be at least 90% for one of the features
# THIS IS STRANDED
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
       -wa \
       -F 0.9 \
       -f 0.9 \
       -s \
       -e \
       -u \
       -a ../gtfs/flair_filter_transcripts_no_chr_novel_genes.sorted.gtf \
       -b ../gtfs/extended_annotation_grch38_only_novel_genes.sorted.gtf \
       > ../overlap_assessment/just_novel_genes_their_overlap_stranded.gtf

#subtract out the overlap and get just their uniquely novel transcripts
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools subtract \
	-s \
        -a ../gtfs/flair_filter_transcripts_no_chr_novel_genes.sorted.gtf \
        -b ../overlap_assessment/just_novel_genes_their_overlap_stranded.gtf \
        > ../overlap_assessment/their_novel_genes_stranded.gtf


# upon merging, the file becomes a bed file format
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools merge \
	-s \
        -i ../overlap_assessment/their_novel_genes_stranded.gtf \
        > ../overlap_assessment/their_novel_genes_stranded.merged.bed

# This creates an additional column of gap lengths to be graphed later on
python3 01a_create_gap_lengths_column.py ../overlap_assessment/their_novel_genes_stranded.merged.bed

##########################################################################




########## ASSESS WITHOUT STRANDEDNESS, REPORT GLINOS TRANSCRIPTS ###########
 
# interesct the gtfs that have only novel genes. Keep full features from glinos data that overlap
# Overlap must be at least 90% for one of the features
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
       -wa \
       -F 0.9 \
       -f 0.9 \
       -e \
       -u \
       -a ../gtfs/flair_filter_transcripts_no_chr_novel_genes.sorted.gtf \
       -b ../gtfs/extended_annotation_grch38_only_novel_genes.sorted.gtf \
       > ../overlap_assessment/just_novel_genes_their_overlap.gtf

#subtract out the overlap and get just their uniquely novel transcripts
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools subtract \
	-a ../gtfs/flair_filter_transcripts_no_chr_novel_genes.sorted.gtf \
	-b ../overlap_assessment/just_novel_genes_their_overlap.gtf \
	> ../overlap_assessment/their_novel_genes.gtf


# upon merging, the file becomes a bed file format
/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools merge \
	-i ../overlap_assessment/their_novel_genes.gtf \
	> ../overlap_assessment/their_novel_genes.merged.bed

# This creates an additional column of gap lengths to be graphed later on
python3 01a_create_gap_lengths_column.py ../overlap_assessment/their_novel_genes.merged.bed

#########################################################################



#get the values to create the venn diagrams
#TODO
#stranded_venn= 
#not_stranded_venn= 




########### CREATE HISTOGRAMS OF GAP_LENGTHS ############################

# create histograms comparing the gap_lengths between transcripts between our data and glinos data, also create a venn diagram
singularity exec /project/mteb223_uksr/singularity_files/transcriptome_long_reads_GTEx.sif Rscript 01b_graph_overlap_and_gaps.R ../overlap_assessment/our_novel_genes_stranded.merged.txt ../overlap_assessment/their_novel_genes_stranded.merged.txt stranded

singularity exec /project/mteb223_uksr/singularity_files/transcriptome_long_reads_GTEx.sif Rscript 01b_graph_overlap_and_gaps.R ../overlap_assessment/our_novel_genes.merged.txt ../overlap_assessment/their_novel_genes.merged.txt not_stranded

########################################################################



# this was used to look at the overlap of all the entries in the Glinos data, not just the new genes
#/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
#	-wa \
#	-wb \
#	-F 0.9 \
#	-f 0.9 \
#	-e \
#	-a extended_annotation_grch38_only_novel.sorted.gtf \
#	-b flair_filter_transcripts_no_chr.sorted.gtf \
#      	> just_novel.gtf


