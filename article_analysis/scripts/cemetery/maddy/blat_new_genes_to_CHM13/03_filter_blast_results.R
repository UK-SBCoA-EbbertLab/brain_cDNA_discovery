library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
blast_output = args[1]


#blat_results <- read_tsv('output_test_back_to_hg38.blast9', col_names = c('Query_id', 'Subject_id', 'perc_ident', 'alignment_length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score'), comment = '#') %>%

blat_results <- read_tsv(blast_output, col_names = c('Query_id', 'Subject_id', 'perc_ident', 'alignment_length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score'), comment = '#') %>%
	separate(Query_id, c('query_chr', 'query_start', 'query_end'), "[^[:alnum:].]+", remove = FALSE) %>%
	#mutate(total_transcript_size = query_end - query_start) %>%
	arrange(desc(alignment_length), desc(perc_ident)) %>%
	group_by(Query_id)

write_tsv(blat_results %>% slice(1:1), paste0(blast_output, "top_match_for_each.tsv"))
write_tsv(blat_results %>% slice(1:10), paste0(blast_output, "top_10_matchs_for_each.tsv"))






