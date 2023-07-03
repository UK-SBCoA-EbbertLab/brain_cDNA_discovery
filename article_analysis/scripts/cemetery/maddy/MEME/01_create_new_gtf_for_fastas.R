library(tidyverse)

#gtf <- read_tsv("../../../data/maddy/MEME/new_rna_UNFILTERED.gtf", col_names = c("chr", "source", "type", "start", "end", "score", "strand", "frame", "info_col"))
#gtf_new_from_known <- read_tsv("../../../data/maddy/MEME/new_exons_from_new_transcripts_in_known_genes.tsv", col_types = "cccddccciic")
#gtf_new_from_new <- read_tsv("../../../data/maddy/MEME/new_exons_from_new_transcripts_in_new_genes.tsv", col_types = "cccddccciic")
gtf_new_from_mito <- read_tsv("../../../data/maddy/MEME/new_exons_from_new_transcripts_mitochondrial_genes.tsv", col_types = "cccddccciic")
#gtf_known_from_known <- read_tsv("../../../data/maddy/MEME/exons_from_known_spliced_genes.tsv", col_types = "cccddccciic")
#gtf_known_from_known_sampled <- read_tsv("../../../data/maddy/MEME/exons_from_known_spliced_genes.tsv", col_types = "cccddccciic") %>%
#	slice_sample(n = 500000) %>%
#	arrange(chr, start, end)
#write_tsv(gtf_known_from_known_sampled, "../../../data/maddy/MEME/exons_from_known_spliced_genes_sampled.tsv")
#gtf_new_from_mito <- read_tsv("../../../data/maddy/MEME/new_exons_from_new_transcripts_mitochondrial_genes.tsv", col_types = "cccddcccic")

separate_info <- function(tib, five_bed, three_bed, exon_bed, extended_exons) {
	tib <- tib %>%
#		extract(info_col, into= c("gene_id", "transcript_id", "exon"), regex = "gene_id \"([A-Za-z0-9]+)\"; transcript_id \"([A-Za-z0-9]+)\"; exon_number \"([0-9]+)\";") %>%
#		filter(type == "exon") %>%
		mutate(bed_start = start - 1) %>% # bed files are 0-based indexing while gtf files are 1-based (see this wiki page :https://en.wikipedia.org/wiki/BED_(file_format)#:~:text=This%20choice%20is,lines%20are%20used.) 
		mutate(three_prime_splice_start = ifelse(strand == "+", bed_start-12, end-1)) %>%
		mutate(three_prime_splice_end = ifelse(strand == "+", bed_start+1, end+12)) %>%
		mutate(five_prime_splice_start = ifelse(strand == "+", end-3, bed_start-6)) %>%
		mutate(five_prime_splice_end = ifelse(strand == "+", end+6, bed_start+3)) %>%
		add_column(score = 1) %>%
#		select(-c(score, frame)) %>%
		mutate(name = paste0(transcript_id, "_", exon_number))

	# for the 5' splice sites, group by transcript_id then filter out those sites that are not actually splice sites, due to them being on the 
	# first exon or last exon
	five_prime_splices <- tib %>%
		group_by(transcript_id) %>%
	       	filter(exon_number != last_exon) %>%
		ungroup() %>%
		select(chr, five_prime_splice_start, five_prime_splice_end, name, score, strand)
	# for the 3' splice sites, group by transcript_id then filter out those sites that are not actually splice sites, due to them being on the 
	# first or last exon
	three_prime_splices <- tib %>% 
		group_by(transcript_id) %>%
		filter(exon_number != 1) %>%
		ungroup() %>%
		select(chr, three_prime_splice_start, three_prime_splice_end, name, score, strand)
	exons <- tib %>% 
		select(chr, start, end, name, score, strand)

	write_tsv(five_prime_splices, five_bed, col_names=FALSE)
	write_tsv(three_prime_splices, three_bed, col_names=FALSE)
	write_tsv(exons, exon_bed, col_names=FALSE)

}

#gtf <- separate_info(gtf, "../../../data/maddy/MEME/five_prime_splice_sites.bed", "../../../data/maddy/MEME/three_prime_splice_sites.bed", "../../../data/maddy/MEME/exons.bed")
#gtf_new_from_known <- separate_info(gtf_new_from_known, "../../../data/maddy/MEME/five_prime_splice_sites_nfk_12s.bed", "../../../data/maddy/MEME/three_prime_splice_sites_nfk_12s.bed", "../../../data/maddy/MEME/exons_nfk_12s.bed")
#gtf_new_from_new <- separate_info(gtf_new_from_new, "../../../data/maddy/MEME/five_prime_splice_sites_nfn_12s.bed", "../../../data/maddy/MEME/three_prime_splice_sites_nfn_12s.bed", "../../../data/maddy/MEME/exons_nfn_12s.bed")
gtf_new_from_mito <- separate_info(gtf_new_from_mito, "../../../data/maddy/MEME/five_prime_splice_sites_nfm_12s.bed", "../../../data/maddy/MEME/three_prime_splice_sites_nfm_12s.bed", "../../../data/maddy/MEME/exons_nfm_12s.bed")
#gtf_known_from_known <- separate_info(gtf_known_from_known, "../../../data/maddy/MEME/five_prime_splice_sites_kfk_12s.bed", "../../../data/maddy/MEME/three_prime_splice_sites_kfk_12s.bed", "../../../data/maddy/MEME/exons_kfk_12s.bed")

