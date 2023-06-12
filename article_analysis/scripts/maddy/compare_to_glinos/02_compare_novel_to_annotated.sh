#!/bin/bash

a_gtf=$1
b_gtf=$2
output_file_prefix=$3

if [ "$6" = "None" ]; then
	intersect="${output_file_prefix}_${4}.intersect.gtf"
	no_overlap="${output_file_prefix}_${4}.no_overlap.gtf"
	no_gene_anno="${output_file_prefix}_${4}.intersect.no_gene_annotation.gtf"


	# interesct Ensembl data with glinos data
	/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
		-s \
		-wa \
		-u \
		-a ${a_gtf} \
		-b ${b_gtf} \
		> "..\\${output_file_prefix}\\${intersect}"

	# get the entried with no overlap (using the -v param)
	/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
                -s \
                -wa \
		-v \
                -a ${a_gtf} \
                -b ${b_gtf} \
                > "..\\${output_file_prefix}\\${no_overlap}"


else
	intersect="${output_file_prefix}_${4}_${6}.intersect.gtf"
	no_overlap="${output_file_prefix}_${4}_${6}.no_overlap.gtf"
	no_gene_anno="${output_file_prefix}_${4}_${6}.intersect.no_gene_annotation.gtf"

	# interesct Ensembl data with glinos data
	/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
		-s \
		-wa \
	        -f ${6} \
	        -F ${6} \
	        -e \
		-u \
		-a ${a_gtf} \
		-b ${b_gtf} \
		> "..\\${output_file_prefix}\\${intersect}"

	# Get the entries with no overlap (using the -v param)
	/project/mteb223_uksr/sequencing_resources/singularity_tools/bin/bedtools intersect \
                -s \
                -wa \
                -f ${6} \
                -F ${6} \
                -e \
		-v \
                -a ${a_gtf} \
                -b ${b_gtf} \
                > "..\\${output_file_prefix}\\${no_overlap}"

fi

awk -v FS='\t' -v OFS='\t' ${5} "..\\${output_file_prefix}\\${intersect}" > "..\\${output_file_prefix}\\${no_gene_anno}"


