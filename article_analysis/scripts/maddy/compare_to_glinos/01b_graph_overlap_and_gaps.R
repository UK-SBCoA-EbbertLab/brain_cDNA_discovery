library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

# open the files that contain the gap lengths for our data and for the Glinos data
our_gaps = read_tsv(args[1], col_name=c("chr", "start", "end", "gap_length")) %>%
	add_column(data_name = "our data")
their_gaps = read_tsv(args[2], col_name=c("chr", "start", "end", "gap_length")) %>%
	add_column(data_name = "their data")

prefix = args[3]

# combine the gap data into 1 tibble, then log transform it
combined <- bind_rows(our_gaps, their_gaps) %>%
	select(data_name, gap_length) %>%
	mutate(log_transformed = log(gap_length))

# plot the untransformed gap length histogram against eachother
ggplot(combined, aes(x=gap_length, color=data_name)) + geom_histogram(alpha=0.5, position="identity")
fileName = paste0(prefix,  "histogram_no_transform.png")
ggsave(fileName)

# plot the transformed gap length histograms against eachother
ggplot(combined, aes(x=log_transformed, color=data_name)) + geom_histogram(alpha=0.5, position="identity")
#ggplot(combined, aes(x=gap_length, color=data_name)) + geom_histogram(alpha=0.5, position="identity") + scale_x_continuous(trans = 'log2')
fileName = paste0(prefix, "histogram_transformation.png")
ggsave(fileName)

#TODO Still need to incorperate the venn diagram
#venn_numbers = str_split(args[3], "|")


