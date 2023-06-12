library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

annotation_tibble <- read_tsv(args[1]) %>%
	mutate(File_name = trimws(File_name, whitespace = "_condensed_threshold.*tsv")) %>%
	gather("Diff", "number", -File_name)

annotation_tibble$Diff <- factor(annotation_tibble$Diff, levels=c("Diff_1bp", "Diff_5bp", "Diff_10bp", "Diff_20bp", "Diff_50bp", "Diff_100bp", "Diff_200bp"))

#annotation_tibble

ggplot(data=annotation_tibble, aes(x=Diff, y=number, group=File_name)) +
	geom_line(aes(color=File_name)) + geom_point(aes(color=File_name)) + theme(legend.position="bottom")

annotation_tibble <- annotation_tibble %>% 
	filter(!grepl('glinos_vs_ensembl', File_name))

ggplot(data=annotation_tibble, aes(x=Diff, y=number, group=File_name)) +
	geom_line(aes(color=File_name)) + geom_point(aes(color=File_name)) + theme(legend.position="bottom")

annotation_tibble_2 <- read_tsv(args[2]) %>%
        mutate(File_name = trimws(File_name, whitespace = "_condensed_threshold.*tsv")) %>%
        gather("Diff", "number", -File_name)

annotation_tibble_2$Diff <- factor(annotation_tibble_2$Diff, levels=c("Diff_1bp", "Diff_5bp", "Diff_10bp", "Diff_20bp", "Diff_50bp", "Diff_100bp", "Diff_200bp"))

#annotation_tibble_2

ggplot(data=annotation_tibble_2, aes(x=Diff, y=number, group=File_name)) +
        geom_line(aes(color=File_name)) + geom_point(aes(color=File_name)) + theme(legend.position="bottom")


