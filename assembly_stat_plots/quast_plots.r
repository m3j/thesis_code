
library(tidyverse)
library(patchwork)
#this file is in the scripts folder
assembly_data <- read_tsv("assembly_summary_2.tsv")

color_scheme = c("#8fb4c7", "#204c7d")

largest_contig_boxplot <- assembly_data %>% ggplot(aes(x= Sample_type, y=Largest_contig, fill=Sample_type)) + 
  geom_boxplot() + 
  labs(title = "Length of the largest contig \nin each assembly") +
  xlab("") +
  ylab("Largest contig (bp)") +
  theme_bw() +
  theme(text = element_text(size=25))+
  scale_fill_manual(values=color_scheme)
largest_contig_boxplot



N50_boxplot <- assembly_data %>% ggplot(aes(x= Sample_type, y=N50, fill=Sample_type)) + 
  geom_boxplot() +
  labs(title = "N50 value for each \nassembly") +
  xlab("") +
  ylab("N50 value") +
  theme_bw() +
  theme(text = element_text(size=25))+
  scale_fill_manual(values=color_scheme)
N50_boxplot

total_length_boxplot <- assembly_data %>% ggplot(aes(x= Sample_type, y=Total_length, fill=Sample_type)) +
  geom_boxplot() +
  labs(title = "Total length of all contigs \nin each assembly") +
  xlab("") +
  ylab("Total contig length (bp)") +
  theme_bw()+
  theme(text = element_text(size=25))+
  scale_fill_manual(values=color_scheme)
total_length_boxplot            

nr_contigs_boxplot <- assembly_data %>% ggplot(aes(x= Sample_type, y=nr_contigs, fill=Sample_type)) +
  geom_boxplot() +
  labs(title = "Number of contigs in each \nassembly") +
  xlab("") +
  ylab("Number of contigs") +
  theme_bw()+
  theme(text = element_text(size=25))+
  scale_fill_manual(values=color_scheme)
nr_contigs_boxplot                                                      

(N50_boxplot + largest_contig_boxplot) / (total_length_boxplot + nr_contigs_boxplot) + plot_layout(guides = "collect")
ggsave("sampletype_assemblystats.png", width = 16, height = 8)
#ggsave(filename = "largest_contig_boxplot.jpeg", plot = largest_contig_boxplot, device=, width = 16, height = 8)

largest_contig.p <- pairwise.wilcox.test(assembly_data$Largest_contig, assembly_data$Sample_type, p.adjust.method = "BH") 

Total_length.p <- pairwise.wilcox.test(assembly_data$Total_length, assembly_data$Sample_type, p.adjust.method = "BH") 

nr_contigs.p <- pairwise.wilcox.test(assembly_data$nr_contigs, assembly_data$Sample_type, p.adjust.method = "BH") 

N50.p <- pairwise.wilcox.test(assembly_data$N50, assembly_data$Sample_type, p.adjust.method = "BH", paired=FALSE) 

# _______________________________________
# assembly statistics in general without 

assembly_data %>%
  summarise(mean = mean(N50), n = n())
