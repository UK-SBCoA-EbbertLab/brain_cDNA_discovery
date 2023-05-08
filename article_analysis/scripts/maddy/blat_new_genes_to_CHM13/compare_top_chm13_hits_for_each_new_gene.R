library(tidyverse)



top_hit_per_gene <- tibble(filename = list.files(path = "new_genes_by_chr/", pattern = "*_new_transcripts_for_new_genes.gtf.fa.blast9top_match_for_each.tsv")) %>%
	rowwise() %>%
	mutate(hits_by_chr = list(read_tsv(paste0("new_genes_by_chr/", filename), col_types = 'cciicdiiiiiiidi'))) %>%
	ungroup() %>%
	unnest(hits_by_chr) %>%
	mutate(total_length = query_end - query_start) %>%
	mutate(dif_in_length = total_length - alignment_length) %>%
	mutate(query_chr = case_when(startsWith(query_chr, 'K') | startsWith(query_chr, 'G') ~ query_chr, TRUE ~ paste0("chr", query_chr))) %>%
	mutate(on_same_chr = case_when(query_chr == Subject_id ~ TRUE, TRUE ~ FALSE)) %>%
	mutate(type = case_when(startsWith(query_chr, 'K') | startsWith(query_chr, 'G') ~ 'scaffold', TRUE ~ 'chromosome')) %>%
	mutate(perfect_match = case_when(perc_ident == 100 & dif_in_length < 1 ~ TRUE, TRUE ~ FALSE)) %>%
	mutate(match_type = case_when(perc_ident == 100 & dif_in_length < 1 ~ 'match_length_and_perc_ident', perc_ident == 100 ~ '100_perc_ident', dif_in_length < 1 ~ 'match_length', TRUE ~ 'neither'))


top_hit_per_gene$match_type <- factor(top_hit_per_gene$match_type, levels=c('match_length_and_perc_ident', '100_perc_ident', 'match_length', 'neither'))

print(top_hit_per_gene, width = Inf)



#	mutate(total_samp_reads = sum(X3)) %>%
#	mutate(perc_chr_samp_reads = X3 * 100 / total_samp_reads) %>%
#	filter(X3 > 0)


chromos <- top_hit_per_gene %>%
	select(query_chr) %>%
	distinct(query_chr) %>%
	pull(query_chr)

print(length(chromos))

colorLevels <- setNames(c("#F8766D", "#F27D52", "#EA842F", "#E18A00", "#D79000", "#CB9600", "#BE9C00", "#B0A100",
			  "#9FA600", "#8CAB00", "#75AF00", "#58B300", "#24B700", "#00BA38", "#00BC58", "#00BE70",
			  "#00C086", "#00C199", "#00C1AB", "#00C0BC", "#00BECC", "#00BBDA", "#00B7E7", "#00B2F3",
			  "#00ACFC", "#00A5FF", "#619CFF", "#8B93FF", "#A989FF", "#C27FFF", "#D575FE", "#E56DF5",
			  "#F066EA", "#F962DD", "#FE61CE", "#FF62BD", "#FF65AC", "#FF6A99", "#FD7084"), levels(as.factor(chromos)))

ggplot(top_hit_per_gene, aes(x=dif_in_length, y=perc_ident, fill=query_chr, color=query_chr)) +
	geom_point() #+
#	scale_fill_manual(values = colorLevels)

ggplot(top_hit_per_gene, aes(x=alignment_length, y=perc_ident, fill=query_chr, color=query_chr)) +
	geom_point() #+

ggplot(top_hit_per_gene, aes(x=dif_in_length)) +
	geom_histogram()

ggplot(top_hit_per_gene, aes(x=perc_ident)) +
	geom_histogram()

ggplot(top_hit_per_gene, aes(x=on_same_chr, y=perc_ident, color=query_chr)) +
#	geom_point()
	geom_jitter()

ggplot(top_hit_per_gene, aes(x=query_chr, y=perc_ident, color=on_same_chr)) +
	geom_jitter()


ggplot(top_hit_per_gene %>% filter(on_same_chr == FALSE), aes(x=alignment_length, y=perc_ident, fill=query_chr, color=query_chr)) +
	geom_point() #+

ggplot(top_hit_per_gene, aes(x=on_same_chr, fill=type)) +
	geom_bar(position = 'dodge') +
	geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(0.9), vjust = -0.5)

ggplot(top_hit_per_gene, aes(x=on_same_chr, y=perc_ident, fill = type)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter()

ggplot(top_hit_per_gene, aes(x=perfect_match, fill = on_same_chr)) +
	geom_bar(position = 'dodge') +
	geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(0.9), vjust = -0.5)

ggplot(top_hit_per_gene, aes(x=match_type, fill = on_same_chr)) +
	geom_bar(position = 'dodge') +
	geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(0.9), vjust = -0.5)



#write_tsv(idxfiles, "testing.tsv")


#ggplot(idxfiles, aes(filename)) + 
#	geom_bar(aes(y = perc_chr_samp_reads, fill=as.factor(X1)), stat = 'identity', position = 'fill')

#ggplot(idxfiles %>% filter(type == 'sorted'), aes(x=X1, y=perc_chr_samp_reads)) + 
#	geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
#	ggtitle('All Reads from BambuGene290099 Region') +
#	facet_wrap(~sample) + 
#	geom_text(aes(label = X3), vjust = -0.5) + 
#	scale_fill_manual(values = colorLevels) +
#	theme(
#	      axis.title.x=element_blank(),
#	      axis.text.x=element_blank(),
#	      axis.ticks.x=element_blank())
#ggsave("perc_reads_from_region.pdf")

#ggplot(idxfiles %>% filter(type == 'primary_only'), aes(x=X1, y=perc_chr_samp_reads)) + 
#	geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
#	ggtitle('Primary Reads from BambuGene290099 Region') +
#	facet_wrap(~sample) + 
#	geom_text(aes(label = X3), vjust = -0.5) + 
#	scale_fill_manual(values = colorLevels) +
#	theme(
#	      axis.title.x=element_blank(),
#	      axis.text.x=element_blank(),
#	      axis.ticks.x=element_blank())
#ggsave("perc_reads_from_region_primary_only.pdf")

#ggplot(idxfiles %>% filter(type == 'filtered_mapq_primary_alignments'), aes(x=X1, y=perc_chr_samp_reads)) + 
#	geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
#	ggtitle('Primary Reads with MAPQ > 10 from BambuGene290099 Region') +
#	facet_wrap(~sample) + 
#	geom_text(aes(label = X3), vjust = -0.5) + 
#	scale_fill_manual(values = colorLevels) +
#	theme(
#	      axis.title.x=element_blank(),
#	      axis.text.x=element_blank(),
#	      axis.ticks.x=element_blank())
#ggsave("perc_reads_from_region_filtered_mapq_primary_alignments.pdf")


#idxfiles$type <- setNames(idxfiles$type, levels(c('filtered_mapq_primary_alignments', 'primary_only', 'sorted')))

#idxplots <- idxfiles %>%
#	group_by(sample) %>%
#	nest() %>%
#	mutate(plot = map2(data, sample,
#			   ~ ggplot(data = .x, aes(x=X1, y=perc_chr_samp_reads)) +
#				   ggtitle(.y) +
#				   geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
#				   facet_wrap(~type) +
#				   geom_text(aes(label = X3), vjust = -0.5) + 
#				   scale_fill_manual(values = colorLevels) +
#				   theme(
#	      				axis.title.x=element_blank(),
#					axis.text.x=element_blank(),
#					axis.ticks.x=element_blank())))

#print(idxplots$plot)

