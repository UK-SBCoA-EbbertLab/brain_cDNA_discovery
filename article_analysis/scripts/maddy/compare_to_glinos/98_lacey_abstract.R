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

# the gtex CPM for transcripts, filtered to include only new transcripts and removing MT transcripts/genes
gtex_counts <- read_tsv("compare_to_glinos/GTEx_GRCh38_version_107_mapq_10_track_reads_quantification/bambu_quant/CPM_transcript.txt") %>%
	filter(grepl("Bambu", TXNAME)) %>%
	pivot_longer(cols = -c(1:2), names_to = "sample", values_to = "CPM") %>%
	rowwise() %>%
	mutate(bam_file = get_samp_id(sample)) %>%
	filter(!(GENEID %in% mt_gene_ids))

# the CPM from our transcripts, filtered to include only new transcripts and removing MT transcripts/genes
our_counts <- read_tsv("compare_to_glinos/ad_vs_ct_pilot_study_february_2023_GRCh38-107_discovery/bambu_discovery/CPM_transcript.txt") %>%
	filter(grepl("Bambu", TXNAME)) %>%
	filter(!(GENEID %in% mt_gene_ids))

# getting just the medically relevant genes from our data. This is going to help us calculate the numbers below
our_counts_med_rel <- our_counts %>% filter(GENEID %in% med_rel_ids)

# number of new transripts from known genes in our data
n_nfk_t <- length(our_counts %>% 
	filter(grepl("ENSG", GENEID)) %>%
	select(TXNAME) %>%
	distinct(TXNAME) %>%
	pull(TXNAME))

# number of new transripts from new gene bodies in our data
n_nfn_t <- length(our_counts %>%
        filter(grepl("Bambu", GENEID)) %>%
        select(TXNAME) %>%
        distinct(TXNAME) %>%
        pull(TXNAME))

# total number of known genes that we found new transcripts for
n_nfk_g <- length(our_counts %>%
        filter(grepl("ENSG", GENEID)) %>%
        select(GENEID) %>%
        distinct(GENEID) %>%
        pull(GENEID))

# total number of new gene bodies that we found new transcripts for
n_nfn_g <- length(our_counts %>%
        filter(grepl("Bambu", GENEID)) %>%
        select(GENEID) %>%
        distinct(GENEID) %>%
        pull(GENEID))

# number of new transripts from known genes in our data -- just medically relevant genes
mrn_nfk_t <- length(our_counts_med_rel %>% 
	filter(grepl("ENSG", GENEID)) %>%
	select(TXNAME) %>%
	distinct(TXNAME) %>%
	pull(TXNAME))

# number of new transripts from new gene bodies in our data -- just medically relevant genes (not actually applicable)
mrn_nfn_t <- length(our_counts_med_rel %>%
        filter(grepl("Bambu", GENEID)) %>%
        select(TXNAME) %>%
        distinct(TXNAME) %>%
        pull(TXNAME))

# total number of known genes that we found new transcripts for -- just medically relevant genes
mrn_nfk_g <- length(our_counts_med_rel %>%
        filter(grepl("ENSG", GENEID)) %>%
        select(GENEID) %>%
        distinct(GENEID) %>%
        pull(GENEID))

# total number of new gene bodies that we found new transcripts for -- just medically relevant genes ( not actually applicable)
mrn_nfn_g <- length(our_counts_med_rel %>%
        filter(grepl("Bambu", GENEID)) %>%
        select(GENEID) %>%
        distinct(GENEID) %>%
        pull(GENEID))

# write all the above numbers to a file
our_totals <- tibble(stat = 
		     c("n_transcripts_nfk", 
		       "n_genes_nfk", 
		       "n_transcripts_nfn", 
		       "n_genes_nfn", 
		       "mrn_transcripts_nfk", 
		       "mrn_genes_nfk", 
		       "mrn_transcripts_nfn", 
		       "mrn_genes_nfn"), 
		     count = 
			     c(n_nfk_t, n_nfk_g, n_nfn_t, n_nfn_g, mrn_nfk_t, mrn_nfk_g, mrn_nfn_t, mrn_nfn_g))
write_tsv(our_totals, "compare_to_glinos/98_lacey_abstract/our_data_totals_stats.tsv")

# reformatting our counts
our_counts <- our_counts %>%
        pivot_longer(cols = -c(1:2), names_to = "sample", values_to = "CPM") %>%
	add_column(tissue_site_detail = "Our Brain Samples")

# getting the total number of new transcripts in each category
nfn_total <- length(gtex_counts %>% ungroup() %>% select(GENEID) %>% filter(grepl("Bambu", GENEID)) %>% distinct(GENEID) %>% pull(GENEID))
nfk_total <- length(gtex_counts %>% ungroup() %>% select(GENEID) %>% filter(grepl("ENSG", GENEID)) %>% distinct(GENEID) %>% pull(GENEID))

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

# creating a grouped brain using 3 brain tissues from GTEX
tissue_brain <- tissues %>% filter(grepl("Brain", tissue_site_detail)) %>% mutate(tissue_site_detail = 'Grouped Brain')

# combinging the counts with the tissue
gtex_brain <- left_join(tissue_brain, gtex_counts) %>%
        select(-c(sample)) %>%
	group_by(TXNAME,GENEID, tissue_site_detail) %>%
        summarise(median_cpm = median(CPM)) %>%
	mutate(med_cpm_gt_0 = median_cpm > 0) %>%
	mutate(med_cpm_gt_1 = median_cpm > 1)
        
# gtex samples, calculating the median cpm for each new transcript by tissue
gtex_samps <- left_join(tissues, gtex_counts) %>% 
	select(-c(sample)) %>%
	group_by(TXNAME,GENEID, tissue_site_detail) %>%
        summarise(median_cpm = median(CPM)) %>%
	mutate(med_cpm_gt_0 = median_cpm > 0) %>%
	mutate(med_cpm_gt_1 = median_cpm > 1)

# our samples, calculating the median cpm for each new transcript
our_brain <- our_counts %>% 
	select(-c(sample)) %>%
	group_by(TXNAME,GENEID, tissue_site_detail) %>%
        summarise(median_cpm = median(CPM)) %>%
	mutate(med_cpm_gt_0 = median_cpm > 0) %>%
	mutate(med_cpm_gt_1 = median_cpm > 1)


# putting all the counts/tissue types together
gtex <- bind_rows(gtex_samps, gtex_brain) %>%
	bind_rows(our_brain) %>% 
	arrange(GENEID, TXNAME, median_cpm)

# same as above, but just medically relevant genes
med_rel_gtex <- gtex %>%
	filter(GENEID %in% med_rel_ids) %>%
	left_join(med_rel_genes %>% rename(GENEID = gene_id))

# write tables to file 
write_tsv(gtex, "compare_to_glinos/98_lacey_abstract/quantification_transcript_by_tissue_median_cpm.tsv")
write_tsv(med_rel_gtex, "compare_to_glinos/98_lacey_abstract/quantification_transcript_by_tissue_median_cpm_med_rel.tsv")

# plotting the expression for each transcript per tissue type
plot_exp <- gtex %>% 
       filter(grepl("ENSG", GENEID)) %>%
       filter(!(GENEID %in% mt_gene_ids)) %>%
       #filter(grepl("Bambu", GENEID)) %>%
       #group_by(tissue_site_detail) %>%
       rowwise() %>%
       mutate(extreme = ifelse(median_cpm > 500, TXNAME, NA))
       #mutate(extreme = ifelse(median_cpm > 30, TXNAME, NA))
ggplot(plot_exp, aes(x = tissue_site_detail, y = median_cpm, color = tissue_site_detail)) +
	geom_boxplot(outlier.shape = NA) + 
	geom_jitter() + 
	geom_text(aes(label=extreme), nudge_y=5)

# summarizing the gtex data - % new transcripts witnessed for each tissue type
gtex_summary <- gtex %>%
	ungroup() %>%
	mutate(status = case_when(startsWith(GENEID, "E") ~ "nfk", startsWith(GENEID, "B") ~ "nfn")) %>%
	group_by(tissue_site_detail, status) %>%
	summarise(n_med_cpm_gt_0 = sum(med_cpm_gt_0), n_med_cpm_gt_1 = sum(med_cpm_gt_1)) %>%
	mutate(perc_med_cpm_gt_0 = case_when(status == "nfk" ~ n_med_cpm_gt_0/n_nfk_t, status == "nfn" ~ n_med_cpm_gt_0/n_nfn_t)) %>%
	mutate(perc_med_cpm_gt_1 = case_when(status == "nfk" ~ n_med_cpm_gt_1/n_nfk_t, status == "nfn" ~ n_med_cpm_gt_1/n_nfn_t))

med_rel_gtex_summary <- med_rel_gtex %>%
	ungroup() %>%
	mutate(status = case_when(startsWith(GENEID, "E") ~ "nfk", startsWith(GENEID, "B") ~ "nfn")) %>%
	group_by(tissue_site_detail, status) %>%
	summarise(n_med_cpm_gt_0 = sum(med_cpm_gt_0), n_med_cpm_gt_1 = sum(med_cpm_gt_1)) %>%
	mutate(perc_med_cpm_gt_0 = case_when(status == "nfk" ~ n_med_cpm_gt_0/mrn_nfk_t, status == "nfn" ~ n_med_cpm_gt_0/mrn_nfn_t)) %>%
	mutate(perc_med_cpm_gt_1 = case_when(status == "nfk" ~ n_med_cpm_gt_1/mrn_nfk_t, status == "nfn" ~ n_med_cpm_gt_1/mrn_nfn_t))


write_tsv(gtex_summary, "compare_to_glinos/98_lacey_abstract/quantification_by_tissue_median_cpm.tsv")
write_tsv(med_rel_gtex_summary, "compare_to_glinos/98_lacey_abstract/quantification_by_tissue_median_cpm_med_rel.tsv")

ggplot(gtex_summary %>% 
       		filter(status == 'nfk') %>% 
       		pivot_longer(cols = c(perc_med_cpm_gt_0, perc_med_cpm_gt_1), names_to = "perc_threshold", values_to = 'perc_n_transcripts'),
	aes(x = tissue_site_detail, y = perc_n_transcripts, fill=perc_threshold)) +
geom_bar(stat = 'identity', position = 'dodge')

ggplot(gtex_summary %>%
                filter(status == 'nfn') %>%
       		pivot_longer(cols = c(perc_med_cpm_gt_0, perc_med_cpm_gt_1), names_to = "perc_threshold", values_to = 'perc_n_transcripts'),
        aes(x = tissue_site_detail, y = perc_n_transcripts, fill=perc_threshold)) +
geom_bar(stat = 'identity', position = 'dodge')


