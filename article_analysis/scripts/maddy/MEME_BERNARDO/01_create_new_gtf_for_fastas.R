#!/usr/bin/Rscript

library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
library(ggplot2)


gtf_new_from_known_start <- read_tsv("../../../data/bernardo/processed/99.other/create_annotations/meme_annotation/new_exon_START_boundaries_from_new_transcripts_in_known_genes.tsv", col_types = "cccddccciic")
gtf_new_from_new_start <- read_tsv("../../../data/bernardo/processed/99.other/create_annotations/meme_annotation/new_exon_START_boundaries_from_new_transcripts_in_new_genes.tsv", col_types = "cccddccciic")
gtf_new_from_mito_start <- read_tsv("../../../data/bernardo/processed/99.other/create_annotations/meme_annotation/new_exon_START_boundaries_from_new_transcripts_mitochondrial_genes.tsv", col_types = "cccddccciic")
gtf_known_from_known_start <- read_tsv("../../../data/bernardo/processed/99.other/create_annotations/meme_annotation/exon_START_boundaries_from_known_spliced_genes.tsv", col_types = "cccddccciic")


gtf_new_from_known_end <- read_tsv("../../../data/bernardo/processed/99.other/create_annotations/meme_annotation/new_exon_END_boundaries_from_new_transcripts_in_known_genes.tsv", col_types = "cccddccciic")
gtf_new_from_new_end <- read_tsv("../../../data/bernardo/processed/99.other/create_annotations/meme_annotation/new_exon_END_boundaries_from_new_transcripts_in_new_genes.tsv", col_types = "cccddccciic")
gtf_new_from_mito_end <- read_tsv("../../../data/bernardo/processed/99.other/create_annotations/meme_annotation/new_exon_END_boundaries_from_new_transcripts_mitochondrial_genes.tsv", col_types = "cccddccciic")
gtf_known_from_known_end <- read_tsv("../../../data/bernardo/processed/99.other/create_annotations/meme_annotation/exon_END_boundaries_from_known_spliced_genes.tsv", col_types = "cccddccciic")


separate_info <- function(tib_start, tib_end, five_bed, three_bed) {
	
        tib_start_pos <- tib_start %>%
		filter(strand == "+") %>%
		mutate(bed_start = start - 1) %>% # bed files are 0-based indexing while gtf files are 1-based (see this wiki page :https://en.wikipedia.org/wiki/BED_(file_format)#:~:text=This%20choice%20is,lines%20are%20used.) 
		mutate(three_prime_splice_start = bed_start-12) %>%
		mutate(three_prime_splice_end =  bed_start+1) %>%
		add_column(score = 1) %>%
		mutate(name = paste0(transcript_id, "_", exon_number))


         tib_start_neg <- tib_start %>%
                 filter(strand == "-") %>%
                 mutate(bed_start = start - 1) %>% # bed files are 0-based indexing while gtf files are 1-based (see this wiki page :https://en.wikipedia.org/wiki/BED_(file_format)#:~:te
                 mutate(five_prime_splice_start = bed_start-6) %>%
                 mutate(five_prime_splice_end = bed_start+3) %>%
                 add_column(score = 1) %>%
                 mutate(name = paste0(transcript_id, "_", exon_number))


	tib_end_pos <- tib_end %>%
		filter(strand == "+") %>%
		mutate(bed_start = start - 1) %>% # bed files are 0-based indexing while gtf files are 1-based (see this wiki page :https://en.wikipedia.org/wiki/BED_(file_format)#:~:text=This%20choice%20is,lines%20are%20used.) 
                mutate(five_prime_splice_start = end-3) %>%
                mutate(five_prime_splice_end = end+6) %>%
		add_column(score = 1) %>%
		mutate(name = paste0(transcript_id, "_", exon_number))


	tib_end_neg <- tib_end %>%
		filter(strand == "-") %>%
		mutate(bed_start = start - 1) %>% # bed files are 0-based indexing while gtf files are 1-based (see this wiki page :https://en.wikipedia.org/wiki/BED_(file_format)#:~:text=This%20choice%20is,lines%20are%20used.) 
		mutate(three_prime_splice_start = end-1) %>%
		mutate(three_prime_splice_end = end+12) %>%
		add_column(score = 1) %>%
		mutate(name = paste0(transcript_id, "_", exon_number))


        ## Concatenate three prime ends
        three_prime_splices = rbind(tib_start_pos, tib_end_neg)

	# for the 3' splice sites, group by transcript_id then filter out those sites that are not actually splice sites, due to them being on the 
	# first or last exon
	three_prime_splices <- three_prime_splices %>% 
		group_by(transcript_id) %>%
		filter(exon_number != 1) %>%
		ungroup() %>%
		select(chr, three_prime_splice_start, three_prime_splice_end, name, score, strand)

	write_tsv(three_prime_splices, three_bed, col_names=FALSE)


        ## Concatenate three prime ends
        five_prime_splices = rbind(tib_start_neg, tib_end_pos)

        # for the 5' splice sites, group by transcript_id then filter out those sites that are not actually splice sites, due to them being on the
        # first or last exon
        five_prime_splices <- five_prime_splices %>%
                group_by(transcript_id) %>%
                filter(exon_number != last_exon) %>%
                ungroup() %>%
                select(chr, five_prime_splice_start, five_prime_splice_end, name, score, strand)

        write_tsv(five_prime_splices, five_bed, col_names=FALSE)
}

gtf_new_from_known <- separate_info(gtf_new_from_known_start, gtf_new_from_known_end, "../../../data/bernardo/processed/MEME/five_prime_splice_sites_nfk_12s.bed", "../../../data/bernardo/processed/MEME/three_prime_splice_sites_nfk_12s.bed")
gtf_new_from_new <- separate_info(gtf_new_from_new_start, gtf_new_from_new_end, "../../../data/bernardo/processed/MEME/five_prime_splice_sites_nfn_12s.bed", "../../../data/bernardo/processed/MEME/three_prime_splice_sites_nfn_12s.bed")
gtf_new_from_mito <- separate_info(gtf_new_from_mito_start, gtf_new_from_mito_end, "../../../data/bernardo/processed/MEME/five_prime_splice_sites_nfm_12s.bed", "../../../data/bernardo/processed/MEME/three_prime_splice_sites_nfm_12s.bed")
gtf_known_from_known <- separate_info(gtf_known_from_known_start, gtf_known_from_known_end, "../../../data/bernardo/processed/MEME/five_prime_splice_sites_kfk_12s.bed", "../../../data/bernardo/processed/MEME/three_prime_splice_sites_kfk_12s.bed")

