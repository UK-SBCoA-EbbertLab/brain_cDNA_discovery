#!/bin/bash

# initialize variables
annotated_gtfs=("../KNOWN_annotation_for_cell_line_mass_spec.gtf")
new_gtfs=("../NEW_annotation_for_cell_line_mass_spec.gtf")
ann_fasta=("../intermediate_files/KNOWN_annotation_for_cell_line_mass_spec.fa")
ann_fasta_no_seq_unavail=("../intermediate_files/KNOWN_annotation_for_cell_line_mass_spec_peptides_no_seq_unavail.fa")
#ann_seq_unavail=("../intermediate_files/KNOWN_annotation_for_cell_line_mass_spec_seq_unavailable.gtf")
new_transcript_fasta=("../intermediate_files/NEW_annotation_for_cell_line_mass_spec_transcripts.fa")
#seq_unavail_fasta=("../intermediate_files/_unavailable_transcripts.fa")
#combined_fasta=("../intermediate_files/combined_MEDIAN_transcripts.fa")
combined_peptide_fasta=("../intermediate_files/NEW_annotation_for_cell_line_mass_spec_3rf.fa")
peptide_database_fasta_tmp=("../intermediate_files/ALL_annotation_for_cell_line_mass_spec.fa.tmp")
peptide_database_fasta=("../ALL_annotation_for_cell_line_mass_spec.fa")


# https://dev.to/meleu/how-to-join-array-elements-in-a-bash-script-303a
joinByChar() {
	local IFS="$1"
	shift
	echo "$*"
}

for i in {0..0};
do

	# grab the gene_id for all the new transcripts.	
	awk -F "\t" '$3 == "transcript" { print $9 }' ${new_gtfs[i]} | \
		awk -F "; " '{ print $1"\t"$2 }' | \
		sed 's/\(gene_id "\)\([0-9a-zA-Z]*\)\("\)\t\(transcript_id "\)\([0-9a-zA-Z]*\)\(";\)/\2\t\5/' > gene_ids_for_new_transcripts.tmp

	gene_ids=($(awk -F '\t' '$1 ~ /ENSG/ {print $1}' gene_ids_for_new_transcripts.tmp)) 
	match_pat=$(joinByChar '|' "${gene_ids[@]}")
#	split -l 500 gene_ids_for_new_transcripts.tmp gene_ids_for_new_transcripts_split


	# Grab all the transcript IDs from the annotated gtf
	# Pass them to biomaRt to grab the peptide sequences
	# save it all to ${ann_fasta}
	awk -F "\t" '$3 == "transcript" { print $9 }' ${annotated_gtfs[i]} | \
		awk -F "; " '{ print $1"\t"$2 }' | \
		sed 's/\(gene_id "\)\([A-Z0-9]*\)\("\)\t\(transcript_id "\)\([0-9A-Z]*\)\(";\)/\2\t\5/' > known_ids_genes.tmp

	# 
	awk -F "\t" '{print $2}' known_ids_genes.tmp > known_transcript_ids_genes.tmp

	# combine the new transcript IDs with the known ones as well
	cat known_ids_genes.tmp gene_ids_for_new_transcripts.tmp > combined_gene_ids_for_headers.tmp

#
	rm -f ${ann_fasta[i]}

	# grab the peptide sequence from biomaRt for the known proteins
	cat known_transcript_ids_genes.tmp | singularity exec /share/singularity/images/ccs/users/mlpa241/wflow-1.sinf Rscript create_biomaRt_annotated_database.R \
		>> ${ann_fasta[i]}

	# remove the sequence unavailable lines from the annotated peptide fasta
	grep -B 1 "^[A-Z\*][A-Z\*]*$" ${ann_fasta[i]} | grep -v "^--$" > ${ann_fasta_no_seq_unavail[i]}

	# create dna fastas for the new gtfs
	singularity exec /share/singularity/images/ccs/users/mlpa241/wflow-1.sinf gffread -F -w ${new_transcript_fasta[i]} -g ../../Homo_sapiens.GRCh38_ERCC.fa ${new_gtfs[i]}

	# generate the peptide fasta
	singularity exec /share/singularity/images/ccs/users/mlpa241/wflow-1.sinf pypgatk_cli.py dnaseq-to-proteindb \
		--config_file  ../pypgatk_ensembl_config.yaml \
		--input_fasta ${new_transcript_fasta[i]} \
		--output_proteindb ${combined_peptide_fasta[i]} \
		--biotype_str '' \
		--num_orfs 3

	#combine the peptide fastas
	cat ${combined_peptide_fasta[i]} ${ann_fasta_no_seq_unavail[i]} > ${peptide_database_fasta_tmp[i]}

	# create the uniprot style headers
	python create_uniprot_headers.py combined_gene_ids_for_headers.tmp "${peptide_database_fasta_tmp[i]}" "${peptide_database_fasta[i]}"

	# add proteins for dr. strom

	echo ">db|circMAN2A1|ENSG00000112893 MAN2A1 OS=Homo sapien OX=9606 GN=MAN2A1" >> "${peptide_database_fasta[i]}"
	echo "MVEFGSKDLTLLMNLMNGTLNPFKSLWCLIPITTQGQLSMLQEKIDHLERLLAENNEIISNIRDSVINLSESVEDGPKSSQSNFSQGAGSHLLPSQLSLSVDTADCLFASQSGSHNSDVQMLDVYSLISFDNPDGGVWKQGFDITYESNEWDTEPLQVFVVPHSHNDPGPALNVARKNRPFGAFAS" >> "${peptide_database_fasta[i]}"
	echo ">db|new_tau_orf|ENSG00000186868 MAPT OS=Homo sapien OX=9606 GN=MAPT" >> "${peptide_database_fasta[i]}"
	echo "MVKRRSPHRGEQPLQARRARPTPPGFQQKPRPLQRHHPALVNLQNQGIAAATAAPAPQALPAAAPAPRPFQPHPPGSPRRWQWSVLHPSRRLPPRAACRQPPCPCQT" >> "${peptide_database_fasta[i]}"



	


done



