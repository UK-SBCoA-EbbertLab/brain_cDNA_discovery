library(tidyverse)

gtf <- read_tsv("new_transcripts_for_new_genes.gtf", col_names = c("seqname","source","feature","start","end","score","strand","frame","attribute")) %>% 
	separate(attribute, c("gene_id", "transcript_id", "exon_number"), sep = ";", remove = TRUE) %>%
        extract(gene_id, c("gene_id"), regex = "gene_id \"(.*)\"") %>%
        extract(transcript_id, c("transcript_id"), regex = " transcript_id \"(.*)\"") %>%
        extract(exon_number, c("exon_number"), regex = " exon_number \"(.*)\"", convert=TRUE) %>%
	select(seqname, start, end, gene_id, transcript_id)

write_tsv(gtf, 'new_transcripts_for_new_genes.key')

