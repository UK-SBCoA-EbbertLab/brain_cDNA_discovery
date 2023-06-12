library(tidyverse)
library(broom)

gtex_counts <- read.table("compare_to_glinos/GTEx_GRCh38_version_107_mapq_10_track_reads_quantification/bambu_quant/CPM_transcript.txt", sep = '\t', header=TRUE)
sample_metadata <- read.table("compare_to_glinos/GTEx_v9_ONT_metadata.txt", sep = '\t', header=TRUE)

get_samp_id <- function(str) {
        name_array <- unlist(strsplit(str, '.', fixed = TRUE))
        last_element <- unlist(strsplit(name_array[length(name_array)], '_', fixed = TRUE))[1]
        string_1 <- paste0(head(name_array, -1), collapse='-')
        string_2 <- paste(last_element, 'bam', sep = ".")
        return(paste(string_1, string_2, sep = "."))
}
gtex_counts_t <- as.data.frame(t(gtex_counts[, -c(1:2)]))
names(gtex_counts_t) <- gtex_counts$TXNAME
rownames(gtex_counts_t) <- names(gtex_counts)[-c(1:2)]
rownames(gtex_counts_t) <- lapply(rownames(gtex_counts_t), get_samp_id)

## TISSUE DATA FOR THE GTEx
tissues <- as_tibble(sample_metadata) %>%
        select(c(3,4)) %>%
	arrange(bam_file) %>%
	mutate(bam_file = as.factor(bam_file))
print(tissues)

gtex_tib <- as_tibble(gtex_counts_t) %>%
	add_column(rownames(gtex_counts_t), .before=1) %>%
	rename(bam_file=`rownames(gtex_counts_t)`) %>%
	arrange(bam_file) %>%
	mutate(bam_file=as.factor(bam_file))

print(gtex_tib)
tis <- as.data.frame(tissues)
gtib <- as.data.frame(gtex_tib)

gtex <- as_tibble(merge(gtib, tis)) %>% 
	mutate(tissue_site_detail = as_factor(tissue_site_detail)) %>%
	filter(tissue_site_detail %in% c('Muscle - Skeletal',
					 'Lung',
					 'Liver',
					 'Heart - Left Ventricle', 'Heart - Atrial Appendage', 'Cells - Cultured fibroblasts')) %>%
#					 'Brain - Putamen (basal ganglia)', 'Brain - Frontal Cortex (BA9)', 'Brain - Cerebellar Hemisphere')) %>%
	pivot_longer(cols = starts_with("Bambu"), names_to = "transcript", values_to = "CPM") %>%
	group_by(transcript) %>%
	nest() %>% 
	mutate(kruskal_test = map(data, function(df) tidy(kruskal.test(CPM ~ tissue_site_detail, data=df)))) %>%
	mutate(pairwise_wilcox_test = map(data, function(df) tidy(pairwise.wilcox.test(df$CPM, df$tissue_site_detail, p.adjust.method = "BH")))) %>%
	unnest(cols = c(kruskal_test, pairwise_wilcox_test), names_repair = 'universal') %>%
	rename(kw_chi_squared = statistic) %>% 
	rename(kw_pvalue=p.value...4) %>%
	rename(kw_degrees_of_freedom = parameter) %>%
	rename(wilcox_pvalue = p.value...9)



testing<-tidy(pairwise.wilcox.test(gtex$BambuTx101, gtex$tissue_site_detail, p.adjust.method = "BH")) %>% pivot_wider(names_from = group2, values_from - p.value) 


#gtex <- merge(gtib, tis)



summary(gtex)
