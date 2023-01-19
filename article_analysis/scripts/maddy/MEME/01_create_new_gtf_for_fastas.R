library(tidyverse)

gtf <- read_tsv("../../../data/maddy/MEME/new_rna_UNFILTERED.gtf", col_names = c("chr", "source", "type", "start", "end", "score", "strand", "frame", "info_col"))

separate_info <- function(tib) {
	tib <- tib %>%
		extract(info_col, into= c("gene_id", "transcript_id", "exon"), regex = "gene_id \"([A-Za-z0-9]+)\"; transcript_id \"([A-Za-z0-9]+)\"; exon_number \"([0-9]+)\";") %>%
		filter(type == "exon") %>%
		mutate(three_prime_splice_start = ifelse(strand == "+", start-12, end)) %>% #end - 1
		mutate(three_prime_splice_end = ifelse(strand == "+", start+1, end+13)) %>%
		mutate(five_prime_splice_start = ifelse(strand == "+", end-2, start-6)) %>%
		mutate(five_prime_splice_end = ifelse(strand == "+", end+7, start+3)) %>%
		select(-c(score, frame)) %>%
		mutate(name = paste0(transcript_id, "_", exon))

}

gtf <- separate_info(gtf)



five_prime_splices <- gtf %>% select(chr, five_prime_splice_start, five_prime_splice_end, name)

three_prime_splices <- gtf %>% select(chr, three_prime_splice_start, three_prime_splice_end, name)

exons <- gtf %>% select(chr, start, end, name)

write_tsv(five_prime_splices, "../../../data/maddy/MEME/five_prime_splice_sites.bed", col_names=FALSE)
write_tsv(three_prime_splices, "../../../data/maddy/MEME/three_prime_splice_sites.bed", col_names=FALSE)
write_tsv(exons, "../../../data/maddy/MEME/exons.bed", col_names=FALSE)
