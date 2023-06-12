library(tidyverse)

# get a list of all the ids that come from the MT -- going to filter these out
mt_gene_ids <- read_tsv("compare_to_glinos/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/bambu_discovery/extended_annotations.gtf", col_names = FALSE) %>%
        filter(X1 == "MT") %>%
        rowwise() %>%
        mutate(GENEID = strsplit(X9, split='"')[[1]][2]) %>%
        pull(GENEID)
mt_gene_ids <- unique(mt_gene_ids)

# load in the list of medically relevant genes and grab the ids to use for filtering later
med_rel_genes <- read_tsv("compare_to_glinos/medically_relevant_genes_02-04-2023_UPDATED.tsv")
med_rel_ids <- med_rel_genes %>% pull(gene_id)


# a funtion to take the sample name and turn it into the bam file name - this helps link up the gtex tissue types and the counts
get_samp_id <- function(str) {
        name_array <- unlist(strsplit(str, '.', fixed = TRUE))
        last_element <- unlist(strsplit(name_array[length(name_array)], '_', fixed = TRUE))[1]
        string_1 <- paste0(head(name_array, -1), collapse='-')
        string_2 <- paste(last_element, 'bam', sep = ".")
        return(paste(string_1, string_2, sep = "."))
}

# the tissue types - filtered
tissues <- read_tsv("compare_to_glinos/GTEx_v9_ONT_metadata.txt") %>%
        select(c(3,4)) %>%
        filter(!grepl("direct", bam_file)) %>%
        arrange(bam_file) %>%
        filter(tissue_site_detail %in% c('Muscle - Skeletal',
                                         'Lung',
                                         'Liver',
                                         'Heart - Left Ventricle', 'Heart - Atrial Appendage', 'Cells - Cultured fibroblasts',
                                         'Brain - Putamen (basal ganglia)', 'Brain - Frontal Cortex (BA9)', 'Brain - Cerebellar Hemisphere'))


full_counts <- read_tsv('compare_to_glinos/GTEx_GRCh38_version_107_mapq_10_track_reads_quantification/bambu_quant/fullLengthCounts_transcript.txt') %>%
	filter(grepl("Bambu", TXNAME)) %>%
	filter(!(GENEID %in% mt_gene_ids)) %>%
	pivot_longer(cols = -c(1:2), names_to = "sample", values_to = "full_counts") %>%
        rowwise() %>%
        mutate(bam_file = get_samp_id(sample)) %>%
	right_join(tissues) %>%
	select(-sample) %>% 
	ungroup() %>%
	group_by(TXNAME, GENEID, tissue_site_detail) %>%
	summarise(total_full_counts = sum(full_counts))

write_tsv(full_counts, "GTEx_summary_full_length_counts_by_transcript_and_tissue_type.tsv")
#	mutate(total = select(., contains("GTEX")) %>% rowSums()) %>%
#	select(TXNAME, GENEID, total)

full_counts <- full_counts %>%
	ungroup() %>%
	filter(total_full_counts >= 10)

summary_normal <- full_counts %>%
        mutate(type = ifelse(startsWith(GENEID, "E"), "nfk", "nfn")) %>%
	group_by(type, tissue_site_detail) %>%
	summarise(n_transcripts = n())


summary_med_rel <- full_counts %>%
	filter(GENEID %in% med_rel_ids) %>%
        mutate(type = ifelse(startsWith(GENEID, "E"), "mrnfk", "mrnfn")) %>%
	group_by(type, tissue_site_detail) %>%
	summarise(n_transcripts = n())

summary <- bind_rows(summary_normal, summary_med_rel)
write_tsv(summary, "fullLengthCounts_stats.tsv")
