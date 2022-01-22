rm(list=ls())

library(tidyverse)
library(broom)
library(phyloseq)
library(DAtest)
library(vegan)
library(stringr)
library(ggthemes)
library(pairwiseAdonis)
library(patchwork)

# getwd()
# setwd("~/")

otutab <- read.csv("chick_kaiju.otu.tab", sep = "\t", row.names = 1, header = TRUE)

otutab_oknames <- otutab %>% 
  rename_at(vars(ends_with("_out")), funs(str_replace(., "_out", "")))

colnames(otutab_oknames)

OTU = otu_table(otutab_oknames, taxa_are_rows = TRUE)

taxtab <- read.csv("chick_kaiju.tax.tab", sep = "\t", row.names = 1, header = TRUE)
taxmat = as.matrix(taxtab)
TAX = tax_table(taxmat)

setdiff(rownames (otutab_oknames), rownames(taxmat))

metadata = read.csv("../Metadata_56samples.tsv", sep = "\t", header=TRUE)

metadata <- metadata %>%  
  mutate ( region = if_else(grepl("china",farm_ID), 'China', 'Europe'))


META = sample_data(metadata)
rownames(META) <-metadata$ID
physeq = phyloseq(OTU, TAX, META)
physeq 

saveRDS(physeq, "chicken.phyloseq.rds")


physeq_df <- psmelt(physeq)
rawkaijubarplot <- ggplot(physeq_df, aes(x = Sample, y = Abundance, fill = Domain)) + 
  theme_bw() +   
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle=90, size=6))
rawkaijubarplot
# dir.create("Results")
# ggsave("Results/rawkaijubarplot.png", width = unit(15,"cm")) 



# Relative Abundance on genus level
cut_off_min.reads = sum(sample_sums(physeq))*0.00005 
physeq_genus_cutoff = preDA(physeq, min.reads = cut_off_min.reads)
glom_or <- tax_glom(physeq_genus_cutoff, taxrank ="Genus")
glom_or

gpt <- subset_taxa(glom_or)
gpt <- prune_taxa(names(sort(taxa_sums(glom_or),TRUE)[1:20]), glom_or)
top25 <- names(sort(taxa_sums(gpt), decreasing=TRUE))[1:20]
ps.top25 <- transform_sample_counts(gpt, function(OTU) 100* OTU/sum(OTU))


vf2=c("peachpuff2","burlywood2","Peachpuff","blanchedalmond", 
      "palegoldenrod","lemonchiffon1", "lightgoldenrod1", "brown1",
      "darkseagreen1",  "seagreen1","darkseagreen3","steelblue3",  
      "steelblue1","brown3","red2","steelblue3",  "steelblue1", "skyblue1",
      "lightsteelblue1", "lightblue2","lightsteelblue", "lightcyan1","coral",
      "sienna1","tan1","orange","brown3","lightsalmon","red2")
col = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
        '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', 
        '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
        '#000075', '#808080', '#ffffff', '#000000', "red", "pink")


pspe <- plot_bar(ps.top25, fill = "Genus") + 
  ylab("Relative abundance %") + 
  geom_bar(stat= "identity", color="gray") + 
  scale_fill_manual(values = vf2) 
# +  
  # theme(axis.ticks.x = element_blank()) + scale_x_discrete(limits = positions)
pspe + theme(panel.grid.major = element_blank(), 
             axis.title.y = element_text( size=21),
             axis.title.x = element_blank(),
             axis.text.x  = element_text(face = "italic", size=16), 
             axis.text.y = element_text(size=22), 
             legend.title = element_text(size=20),
             legend.text = element_text(face = "italic", size=22), 
             legend.key.size = unit(0.5, "cm"),
             panel.grid.minor = element_blank(), 
             panel.background = element_blank())

ggsave("Results/genus_abundance20.png", width = unit(15,"cm")) 


# Relative Abundance on Family level
physeq_family_cutoff = preDA(physeq, min.reads = cut_off_min.reads)
glom_or_fam <- tax_glom(physeq_family_cutoff, taxrank ="Family")


gpt_fam <- subset_taxa(glom_or_fam)
gpt_fam <- prune_taxa(names(sort(taxa_sums(glom_or_fam),TRUE)[1:15]), glom_or_fam)
top15_fam <- names(sort(taxa_sums(gpt_fam), decreasing=TRUE))[1:15]
ps.top15_fam <- transform_sample_counts(gpt_fam, function(OTU) 100* OTU/sum(OTU))

col15 = c("burlywood2","Peachpuff","blanchedalmond", 
          "palegoldenrod","lemonchiffon1", "brown1",
          "darkseagreen1","darkseagreen3",
          "steelblue1","brown3","lightsalmon","steelblue3",
          "lightsteelblue1","coral",
          "sienna1")

pspe_fam_15 <- plot_bar(ps.top15_fam, fill = "Family") + ylab("Relative abundance %") + 
  geom_bar(stat= "identity", color="gray") + scale_fill_manual(values = col15) 
# +  
# theme(axis.ticks.x = element_blank()) + scale_x_discrete(limits = positions)
pspe_fam_15 + theme(panel.grid.major = element_blank(), 
             axis.title.y = element_text( size=21),
             axis.title.x = element_blank(),
             axis.text.x  = element_text(face = "italic", size=16), 
             axis.text.y = element_text(size=22), 
             legend.title = element_text(size=20),
             legend.text = element_text(face = "italic", size=22), 
             legend.key.size = unit(0.5, "cm"),
             panel.grid.minor = element_blank(), 
             panel.background = element_blank())

ggsave("Results/family_abundance_top15.png", width = unit(15,"cm")) 






# PCA - on all samples and OTUs from bacteria
PCA_color_scheme = c("#a6cee3","#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                     "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", 
                     "#fff099", "#b15928")

physeq_bac <- subset_taxa(physeq, Domain == "Bacteria")
# filter taxa and keep those with at least 0.00005 fraction of the total reads:
cut_off_min.reads = sum(sample_sums(physeq_bac))*0.00005 
physeq_bac_cutoff = preDA(physeq_bac, min.reads = cut_off_min.reads)
physeq_bac_cutoff

#CRL transformation
# zero correction by substituting 0 with 1 and all non-zero values 
# are corrected such that log-ratios before and after correction are the same.
physeq_zc <- transform_sample_counts(physeq_bac_cutoff, 
                                     function(y) sapply(y, function(x) 
                                       ifelse(x==0, 1, (1-(sum(y==0)*1)/sum(y))*x)))
# log-transform
physeq_clr <- transform_sample_counts(physeq_zc, function(x) log(x/exp(mean(log(x))))) 

# redundance analysis RDA without restrictions are the same as PCA.
ord_clr <- ordinate(physeq_clr, "RDA")

plot_scree(ord_clr) +
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
# ggsave("Results/clrscree.png")

clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)

all_PCA_plot <- plot_ordination(physeq, ord_clr, title = "PCA on full dataset", 
                                type="samples", color="farm_ID", shape = "region") +
  geom_point(size = 4) +
  scale_color_manual(values=PCA_color_scheme) + 
  coord_fixed(clr2 / clr1) +
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1))
  # geom_text(aes(label=ID), colour="black")
# ggsave("Results/clrPC_low_without_label.png", width = unit(15, "cm"))



# PCA only related to lactobacillaceae species: 
physeq_Lactobacillaceae<- subset_taxa(physeq, Family == " Lactobacillaceae")

cut_off_min.reads = sum(sample_sums(physeq_Lactobacillaceae))*0.00005 
physeq_Lacto_cutoff = preDA(physeq_Lactobacillaceae, min.reads = cut_off_min.reads)
physeq_Lacto_cutoff

#CRL transformation
# zero correction by substituting 0 with 1 and all non-zero values 
# are corrected such that log-ratios before and after correction are the same.
physeq_Lacto_zc<- transform_sample_counts(physeq_Lacto_cutoff, 
                                     function(y) sapply(y, function(x) 
                                       ifelse(x==0, 1, (1-(sum(y==0)*1)/sum(y))*x)))
# log-transform
physeq_lacto_clr <- transform_sample_counts(physeq_Lacto_zc, function(x) log(x/exp(mean(log(x))))) 
lacto_ord_clr <- ordinate(physeq_lacto_clr, "RDA")

plot_scree(lacto_ord_clr) +
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
# ggsave("Results/clrscree_Lactobacillaceae.png")

lacto_clr1 <- lacto_ord_clr$CA$eig[1] / sum(lacto_ord_clr$CA$eig)
lacto_clr2 <- lacto_ord_clr$CA$eig[2] / sum(lacto_ord_clr$CA$eig)

lacto_PCA_plot <- plot_ordination(physeq, lacto_ord_clr, 
                                  title = "PCA on Lactobacillaceae", type="samples", 
                                color="farm_ID", shape = "region") +
  geom_point(size = 4) +
  scale_color_manual(values=PCA_color_scheme) +
  coord_fixed(lacto_clr2 / lacto_clr1) +
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1),
        legend.position = 'none') 
  # geom_text(aes(label=ID), colour="black")
# ggsave("Results/clrPC_Lactobacillaceae_low_withoutlabel.png", width = unit(15, "cm"))

all_PCA_plot / lacto_PCA_plot
# ggsave("Results/PCAboth.png", width = unit(15, "cm"))

# significant difference in clusters? All
clr_dist_matrix <- distance(physeq_clr, method = "euclidean")
adonis(clr_dist_matrix ~ sample_data(physeq_clr)$region, method = "eucledian")

# significant difference in clusters? Lactobacillaceae
lacto_clr_dist_matrix <- distance(physeq_lacto_clr, method = "euclidean")
adonis(lacto_clr_dist_matrix ~ sample_data(physeq_lacto_clr)$region, method = "eucledian")


# significant difference in FARM clusters? All
clr_dist_matrix <- distance(physeq_clr, method = "euclidean")
adonis(clr_dist_matrix ~ sample_data(physeq_clr)$farm_ID, method = "eucledian")
pairwise.adonis(clr_dist_matrix, sample_data(physeq_clr)$farm_ID, sim.method = "eucledian",
                p.adjust.m = "BH")

# significant difference in FARM clusters? Lactobacillaceae
lacto_clr_dist_matrix <- distance(physeq_lacto_clr, method = "euclidean")
adonis(lacto_clr_dist_matrix ~ sample_data(physeq_lacto_clr)$farm_ID, method = "eucledian")
pairwise.adonis(lacto_clr_dist_matrix, sample_data(physeq_lacto_clr)$farm_ID,
                sim.method = "eucledian", p.adjust.m = "BH")



# plot of the difference of reads percentage of different OTUs within Lactobacillaceae
# looks like kaiju only finds lactubacillus and weissella so 
# just look at the family level of lactobacillaceae.
physeq_lacto_cutoff_df <- psmelt(physeq_Lacto_cutoff)

totalreads <- colSums(otutab_oknames)

#adding the total number of reads for each sample. so we can get the precentage abundance
physeq_lacto_cutoff_df <- physeq_lacto_cutoff_df %>% mutate(total_reads = case_when(ID == "ERR2241649" ~ totalreads["ERR2241649"],
                                                                ID == "ERR2241919" ~ totalreads["ERR2241919"],
                                                                ID == "ERR2241893" ~ totalreads["ERR2241893"],
                                                                ID == "ERR2241925" ~ totalreads["ERR2241925"], 
                                                                ID == "ERR2241887" ~ totalreads["ERR2241887"],
                                                                ID == "ERR2241848" ~ totalreads["ERR2241848"],
                                                                ID == "ERR2241655" ~ totalreads["ERR2241655"],
                                                                ID == "ERR2241647" ~ totalreads["ERR2241647"],
                                                                ID == "ERR2241860" ~ totalreads["ERR2241860"],
                                                                ID == "ERR2241653" ~ totalreads["ERR2241653"],
                                                                ID == "ERR2241890" ~ totalreads["ERR2241890"],
                                                                ID == "ERR2241923" ~ totalreads["ERR2241923"],
                                                                ID == "ERR2241849" ~ totalreads["ERR2241849"],
                                                                ID == "ERR2241945" ~ totalreads["ERR2241945"],
                                                                ID == "SRR6323384" ~ totalreads["SRR6323384"],
                                                                ID == "SRR6323197" ~ totalreads["SRR6323197"],
                                                                ID == "ERR2241899" ~ totalreads["ERR2241899"],
                                                                ID == "SRR6323414" ~ totalreads["SRR6323414"],
                                                                ID == "SRR6323530" ~ totalreads["SRR6323530"],
                                                                ID == "ERR2241891" ~ totalreads["ERR2241891"],
                                                                ID == "SRR6323338" ~ totalreads["SRR6323338"],
                                                                ID == "SRR9113725" ~ totalreads["SRR9113725"],
                                                                ID == "ERR2241927" ~ totalreads["ERR2241927"],
                                                                ID == "SRR6323350" ~ totalreads["SRR6323350"],
                                                                ID == "SRR6323172" ~ totalreads["SRR6323172"],
                                                                ID == "SRR9113737" ~ totalreads["SRR9113737"],
                                                                ID == "ERR2241855" ~ totalreads["ERR2241855"],
                                                                ID == "SRR6323252" ~ totalreads["SRR6323252"],
                                                                ID == "ERR2241856" ~ totalreads["ERR2241856"],
                                                                ID == "ERR2241651" ~ totalreads["ERR2241651"],
                                                                ID == "SRR6323516" ~ totalreads["SRR6323516"],
                                                                ID == "SRR9113707" ~ totalreads["SRR9113707"],
                                                                ID == "SRR9113691" ~ totalreads["SRR9113691"],
                                                                ID == "SRR6323243" ~ totalreads["SRR6323243"],
                                                                ID == "SRR6323514" ~ totalreads["SRR6323514"],
                                                                ID == "SRR6323503" ~ totalreads["SRR6323503"],
                                                                ID == "SRR6323367" ~ totalreads["SRR6323367"],
                                                                ID == "SRR6323556" ~ totalreads["SRR6323556"],
                                                                ID == "SRR6323296" ~ totalreads["SRR6323296"],
                                                                ID == "SRR6323205" ~ totalreads["SRR6323205"],
                                                                ID == "SRR6323295" ~ totalreads["SRR6323295"],
                                                                ID == "SRR6323550" ~ totalreads["SRR6323550"],
                                                                ID == "SRR6323502" ~ totalreads["SRR6323502"],
                                                                ID == "SRR9113706" ~ totalreads["SRR9113706"],
                                                                ID == "SRR6323244" ~ totalreads["SRR6323244"],
                                                                ID == "SRR6323167" ~ totalreads["SRR6323167"],
                                                                ID == "SRR6323493" ~ totalreads["SRR6323493"],
                                                                ID == "SRR6323200" ~ totalreads["SRR6323200"],
                                                                ID == "SRR6323247" ~ totalreads["SRR6323247"],
                                                                ID == "SRR6323343" ~ totalreads["SRR6323343"],
                                                                ID == "SRR6323088" ~ totalreads["SRR6323088"],
                                                                ID == "SRR6323449" ~ totalreads["SRR6323449"],
                                                                ID == "SRR6323134" ~ totalreads["SRR6323134"],
                                                                ID == "SRR6323340" ~ totalreads["SRR6323340"],
                                                                ID == "SRR6323249" ~ totalreads["SRR6323249"],
                                                                ID == "SRR6323250" ~ totalreads["SRR6323250"]
))
head(physeq_lacto_cutoff_df)

physeq_lacto_cutoff_df <- physeq_lacto_cutoff_df %>% mutate(percent_abundance = Abundance/total_reads*100) 


abun_percent <- ggplot(physeq_lacto_cutoff_df, 
                       mapping = aes(x = percent_abundance,  fill = region)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("darkseagreen","lightblue3")) +
  coord_flip() +
  xlab("Abundance in %") +
  theme_set(theme_bw())  +
  theme(text = element_text(size=20),
        legend.text = element_text(face = "italic", size=22),
        axis.text.x = element_blank())
  # theme(axis.ticks.x=element_blank(), 
  #       axis.text.x=element_blank()) 


physeq_lacto_cutoff_P <- pairwise.wilcox.test(physeq_lacto_cutoff_df$percent_abundance, physeq_lacto_cutoff_df$region, p.adjust.method = "BH") 


physeq_lacto_cutoff_df %>% summarise(mean = mean(percent_abundance), n = n())


physeq_lacto_cutoff_filtered <- physeq_lacto_cutoff_df %>%
  filter(percent_abundance < 0.05)


abun_percent_filtered <- ggplot(physeq_lacto_cutoff_filtered, 
                                mapping = aes(x = percent_abundance,  fill = region)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("darkseagreen","lightblue3")) +
  coord_flip() +
  xlab("Abundance in %") +
  theme_set(theme_bw())  +
  theme(text = element_text(size=20),
        legend.text = element_text(face = "italic", size=22),
        axis.text.x = element_blank())
  


abun_percent + abun_percent_filtered + plot_layout(guides = "collect")

ggsave("Results/abundancePercentage_region.png", width = unit(10, "cm"))

