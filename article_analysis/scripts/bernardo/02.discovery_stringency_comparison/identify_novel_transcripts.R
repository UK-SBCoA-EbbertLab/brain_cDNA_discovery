library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# get the start position for the transcript
get_transcript_start <- function(data_tib) {
	data_tib <- data_tib %>% filter(feature == 'transcript')
	data_tib$start

}

# get the end position for the transcript
get_transcript_end <- function(data_tib) {
	data_tib <- data_tib %>% filter(feature == 'transcript')
	data_tib$end
}


# calculate the bp difference for the transcripts
calculate_diff <- function(overlap_tib, data_ours) {
	# rename their transcript start and end positions
	# fix the order of the exons so that they will be the same between our data and theirs
	# also drop the transcript lines to that we only have the exons
	overlap_tib <- overlap_tib %>%
		rename(t_start = start, t_end = end) %>%
		drop_na(exon_number) %>%
		arrange(t_start) %>%
		rowid_to_column("exon_number_ordered") %>%
		select(t_start, t_end, exon_number_ordered)

	data_ours <- data_ours %>% 
		drop_na(exon_number) %>%
		arrange(start) %>%
		rowid_to_column("exon_number_ordered") %>%
		select(start, end, exon_number_ordered)

	# join our data with theirs on the exon_number
	# then calculate the base pair difference for each exon
	# sum all the base pair differences and return that number
	combined_tib <- full_join(data_ours, overlap_tib) %>%
		rowwise() %>%
		mutate(bp_diff = if_else(start > t_end | end < t_start, abs(start-end), abs(start-t_start) + abs(end-t_end))) %>%
		ungroup() %>%
		summarise(comb_bp_diff = sum(bp_diff, na.rm=TRUE))
	
	return(combined_tib$comb_bp_diff)
}

# comparing transcripts to see which of ours are novel
compare_transcripts <- function(seqname_o, strand_o, transcript_id_o, data_o, exons_o, t_start_o, t_end_o, their_tib, bp_thresh_array) {
	overlap <- their_tib %>% filter(seqname == seqname_o, strand == strand_o) %>%
		filter(t_start <= t_end_o & t_end >= t_start_o)

	# if none of their transcripts even overlap ours, then our transcript is novel
	if (nrow(overlap) == 0) {
		return(paste(rep("NOVEL", length(bp_thresh_array)), collapse=","))
	}

	overlap <- overlap %>% filter(exons == exons_o)

	# if their transcripts don't have the same number of exons, then our transcript is novel
	if (nrow(overlap) == 0) {
		return(paste(rep("NOVEL", length(bp_thresh_array)), collapse=","))
	}

	bp_diff_tib = overlap %>% 
		mutate( bp_diff = map(data, calculate_diff, data_o))

	# if the transcripts have more than a 20 bp difference, then our transcript is novel 
	n = length(bp_thresh_array)
	novel_array = c()
	for (i in 1:n) {
		if (any(bp_diff_tib$bp_diff > as.integer(bp_thresh_array[i]))) {
			novel_array=append(novel_array, "NOVEL")
		} else {
			novel_array=append(novel_array, args[4])
		}
	}

	# if our transcript isn't novel, then it has already been annotated
	return(paste(novel_array, collapse=","))

}

# A function to compare the two tibbles and lable the unique novel transcripts in our data 
# the last two functions pair down the tibble so we only have the transcript rows and not the exon rows, and selects the important columns
compare_data <- function(our_tib, their_tib, bp_thresh_array) {
	thresh_column_names=c()
	for (i in length(bp_thresh_array)) {
		thresh_column_names = append(thresh_column_names, paste(bp_thresh_array, "bp_thresh", sep="_"))
	}

	our_tib %>%
		rowwise() %>%
		mutate(novel = compare_transcripts(seqname, strand, transcript_id, data, exons, t_start, t_end, their_tib, bp_thresh_array)) %>%
		separate(novel, thresh_column_names, sep = ",", remove=TRUE) %>% 
		unnest(cols = c(data)) %>%
		filter(feature == "transcript") %>%
		select(seqname, strand, gene_id, transcript_id, start, end, exons, ends_with("bp_thresh"))

}

# we read in the gtf from our novel transcripts, pull out the transcript_id, gene_id and exon_number into their own columns. Then we nest the tibble by seqname, transcript_id, and strand
# We then create a column that tells us how many exons are in a transcript, the start position of the transcript and the end position of the transcript
# This allows us to first compare the transcripts overall to each other, before we compare transcripts on the exon level.
our_tibble = read_tsv(args[1], col_names = c("seqname","source","feature","start","end","score","strand","frame","attribute")) %>%
	separate(attribute, c("gene_id", "transcript_id", "exon_number"), sep = ";", remove = TRUE) %>%
	extract(gene_id, c("gene_id"), regex = "gene_id \"(.*)\"") %>%
	extract(transcript_id, c("transcript_id"), regex = " transcript_id \"(.*)\"") %>%
	extract(exon_number, c("exon_number"), regex = " exon_number \"(.*)\"", convert=TRUE) %>%
	group_by(seqname, transcript_id, strand) %>%
	nest() %>%
	mutate(
	       exons = map_int(data, nrow) -1,
	       t_start = map_dbl(data, get_transcript_start),
	       t_end = map_dbl(data, get_transcript_end)
	) 

# the glinos data transcripts, formatted the same as our novel transcripts
their_tibble = read_tsv(args[2], col_names =c("seqname","source","feature","start","end","score","strand","frame","attribute")) %>%
	separate(attribute, c("gene_id", "transcript_id", "exon_number"), sep = ";", remove = TRUE) %>%
	extract(gene_id, c("gene_id"), regex = "gene_id \"(.*)\"") %>%
	extract(transcript_id, c("transcript_id"), regex = " transcript_id \"(.*)\"") %>%
	extract(exon_number, c("exon_number"), regex = " exon_number \"(.*)\"", convert=TRUE) %>%
        group_by(seqname, transcript_id, strand) %>%
        nest() %>%
        mutate(
               exons = map_int(data, nrow) -1,
               t_start = map_dbl(data, get_transcript_start),
               t_end = map_dbl(data, get_transcript_end)
        )

our_novel_stuff <- compare_data(our_tibble, their_tibble, c(20, 100, 200))
write_tsv(our_novel_stuff, args[3])
