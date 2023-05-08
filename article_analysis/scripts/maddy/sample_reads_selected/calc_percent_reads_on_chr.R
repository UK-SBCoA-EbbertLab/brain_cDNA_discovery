library(tidyverse)


#sample_579_PAG75663_mapped_filtered_sorted.bam_GL000214_125300_133400_region.fastq.filtered_mapq_primary_alignments.idxstat
#sample_579_PAG75663_mapped_filtered_sorted.bam_GL000214_125300_133400_region.fastq.primary_only.bam


#df <- list.files(path = "../../../data/maddy/sample_reads_selected/", pattern = "*.idxstat") %>%
#	map_df(~read_csv(.) %>% add_column(file = 'test'))
#df

idxfiles <- tibble(filename = list.files(path = "../../../data/maddy/sample_reads_selected/", pattern = "*125300_133400_region.fastq.*.idxstat")) %>%
	rowwise() %>%
	mutate(chrom_reads = list(read_tsv(paste0("../../../data/maddy/sample_reads_selected/", filename), col_names = FALSE))) %>%
	ungroup() %>%
	unnest(chrom_reads) %>%
	extract(filename, c('sample', 'type') , "sample_([0-9]+_[A-Z0-9]+).*fastq.(.*).idxstat") %>%
	group_by(sample, type) %>%
	mutate(total_samp_reads = sum(X3)) %>%
	mutate(perc_chr_samp_reads = X3 * 100 / total_samp_reads) %>%
	filter(X3 > 0)


flagstats <- tibble(filename = list.files(path = "../../../data/maddy/sample_reads_selected/", pattern = "*125300_133400_region.fastq.sorted.flagstat.tsv")) %>%
	rowwise() %>%
	mutate(flag_info = list(read_tsv(paste0("../../../data/maddy/sample_reads_selected/", filename), col_names = FALSE))) %>%
	ungroup() %>%
	unnest(flag_info) %>%
	extract(filename, 'sample', "sample_([0-9]+_[A-Z0-9]+).*") %>%
	filter(X3 == "primary" | X3 == "secondary" | X3 == "supplementary")

flagstats$X1 <- as.numeric(flagstats$X1)

flagstats <- flagstats %>%
	group_by(sample) %>%
	mutate(total_reads = sum(X1)) %>%
	mutate(perc_reads = X1 * 100 / total_reads)
	
chromos <- idxfiles %>%
	select(X1) %>%
	distinct(X1) %>%
	pull(X1)

colorLevels <- setNames(c("#F8766D", "#EC8239", "#DB8E00", "#C79800", 
			  "#AEA200", "#8FAA00", "#64B200", "#00B81B", 
			  "#00BD5C", "#00C085", "#00C1A7", "#00BFC4", 
			  "#00BADE", "#00B2F3", "#00A6FF", "#7C96FF", 
			  "#B385FF", "#D874FD", "#EF67EB", "#FD61D3", 
			  "#FF63B6", "#FF6B94"), levels(as.factor(chromos)))

idxfiles$type <- factor(idxfiles$type, levels=c('sorted', 'primary_only', 'filtered_mapq_primary_alignments'))

write_tsv(idxfiles, "testing.tsv")

print(flagstats)

#ggplot(idxfiles, aes(filename)) + 
#	geom_bar(aes(y = perc_chr_samp_reads, fill=as.factor(X1)), stat = 'identity', position = 'fill')

ggplot(idxfiles %>% filter(type == 'sorted'), aes(x=X1, y=perc_chr_samp_reads)) + 
	geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
	ggtitle('All Reads from BambuGene290099 Region') +
	facet_wrap(~sample) + 
	geom_text(aes(label = X3), vjust = -0.5) + 
	scale_fill_manual(values = colorLevels) +
	theme(
	      axis.title.x=element_blank(),
	      axis.text.x=element_blank(),
	      axis.ticks.x=element_blank())
ggsave("perc_reads_from_region.pdf")

ggplot(idxfiles %>% filter(type == 'primary_only'), aes(x=X1, y=perc_chr_samp_reads)) + 
	geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
	ggtitle('Primary Reads from BambuGene290099 Region') +
	facet_wrap(~sample) + 
	geom_text(aes(label = X3), vjust = -0.5) + 
	scale_fill_manual(values = colorLevels) +
	theme(
	      axis.title.x=element_blank(),
	      axis.text.x=element_blank(),
	      axis.ticks.x=element_blank())
ggsave("perc_reads_from_region_primary_only.pdf")

ggplot(idxfiles %>% filter(type == 'filtered_mapq_primary_alignments'), aes(x=X1, y=perc_chr_samp_reads)) + 
	geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
	ggtitle('Primary Reads with MAPQ > 10 from BambuGene290099 Region') +
	facet_wrap(~sample) + 
	geom_text(aes(label = X3), vjust = -0.5) + 
	scale_fill_manual(values = colorLevels) +
	theme(
	      axis.title.x=element_blank(),
	      axis.text.x=element_blank(),
	      axis.ticks.x=element_blank())
ggsave("perc_reads_from_region_filtered_mapq_primary_alignments.pdf")


#idxfiles$type <- setNames(idxfiles$type, levels(c('filtered_mapq_primary_alignments', 'primary_only', 'sorted')))

idxplots <- idxfiles %>%
	group_by(sample) %>%
	nest() %>%
	mutate(plot = map2(data, sample,
			   ~ ggplot(data = .x, aes(x=X1, y=perc_chr_samp_reads)) +
				   ggtitle(.y) +
				   geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
				   facet_wrap(~type) +
				   geom_text(aes(label = X3), vjust = -0.5) + 
				   scale_fill_manual(values = colorLevels) +
				   theme(
	      				axis.title.x=element_blank(),
					axis.text.x=element_blank(),
					axis.ticks.x=element_blank())))

print(idxplots$plot)

#ggsave(paste0("perc_reads_from_region_", idxplots$sample, ".pdf"), idxplots$plot)

# TOTAL READS
ggplot(idxfiles %>% filter(type == 'sorted'), aes(x=X1, y=X3)) + 
	geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
	ggtitle('All Reads from BambuGene290099 Region') +
	facet_wrap(~sample) + 
	geom_text(aes(label = X3), vjust = -0.5) + 
	scale_fill_manual(values = colorLevels) +
	theme(
	      axis.title.x=element_blank(),
	      axis.text.x=element_blank(),
	      axis.ticks.x=element_blank())
ggsave("reads_from_region.pdf")

ggplot(idxfiles %>% filter(type == 'primary_only'), aes(x=X1, y=X3)) + 
	geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
	ggtitle('Primary Reads from BambuGene290099 Region') +
	facet_wrap(~sample) + 
	geom_text(aes(label = X3), vjust = -0.5) + 
	scale_fill_manual(values = colorLevels) +
	theme(
	      axis.title.x=element_blank(),
	      axis.text.x=element_blank(),
	      axis.ticks.x=element_blank())
ggsave("reads_from_region_primary_only.pdf")

ggplot(idxfiles %>% filter(type == 'filtered_mapq_primary_alignments'), aes(x=X1, y=X3)) + 
	geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
	ggtitle('Primary Reads with MAPQ > 10 from BambuGene290099 Region') +
	facet_wrap(~sample) + 
	geom_text(aes(label = X3), vjust = -0.5) + 
	scale_fill_manual(values = colorLevels) +
	theme(
	      axis.title.x=element_blank(),
	      axis.text.x=element_blank(),
	      axis.ticks.x=element_blank())
ggsave("reads_from_region_filtered_mapq_primary_alignments.pdf")


#idxfiles$type <- setNames(idxfiles$type, levels(c('filtered_mapq_primary_alignments', 'primary_only', 'sorted')))

idxplots <- idxfiles %>%
	group_by(sample) %>%
	nest() %>%
	mutate(plot = map2(data, sample,
			   ~ ggplot(data = .x, aes(x=X1, y=X3)) +
				   ggtitle(.y) +
				   geom_bar(aes(fill = as.factor(X1)), stat = 'identity') +
				   facet_wrap(~type) +
#				   geom_text(aes(label = X3), vjust = -0.5) + 
				   scale_fill_manual(values = colorLevels) +
				   theme(
	      				axis.title.x=element_blank(),
					axis.text.x=element_blank(),
					axis.ticks.x=element_blank())))

print(idxplots$plot)


ggplot(flagstats, aes(sample)) +
	geom_bar(aes(y = perc_reads, fill = as.factor(X3)), stat = 'identity', position = 'fill') + 
	coord_flip()

ggsave("perc_reads_type.pdf")
	

