rm(list=ls())


library(tidyverse)
library(broom)
library(phyloseq)
library(DAtest)
library(vegan)
library(stringr)
library(ggthemes)

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

META = sample_data(metadata)
rownames(META) <-metadata$ID
physeq = phyloseq(OTU, TAX, META)
physeq 

saveRDS(physeq, "chicken.phyloseq.rds")


physeq_df <- psmelt(physeq)
rawkaijubarplot <- ggplot(physeq_df, aes(x = Sample, y = Abundance, fill = Domain)) + theme_bw() +   geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle=90, size=6))
rawkaijubarplot
dir.create("Results")
ggsave("Results/rawkaijubarplot.png", width = unit(15,"cm")) 



# Relative Abundance on genus level

physeq_genus_cutoff = preDA(physeq, min.reads = 1000)
glom_or <- tax_glom(physeq_genus_cutoff, taxrank ="Genus")
glom_or

gpt <- subset_taxa(glom_or)
gpt <- prune_taxa(names(sort(taxa_sums(glom_or),TRUE)[1:30]), glom_or)
top25 <- names(sort(taxa_sums(gpt), decreasing=TRUE))[1:30]
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


pspe <- plot_bar(ps.top25, fill = "Genus") + ylab("Relative abundance %") + 
  geom_bar(stat= "identity", color="gray") + scale_fill_manual(values = vf2) 
# +  
  # theme(axis.ticks.x = element_blank()) + scale_x_discrete(limits = positions)
pspe + theme(panel.grid.major = element_blank(), 
             axis.title.y =element_text(face = "bold", size=11),
             axis.text.x  = element_text(face = "italic", size=12), 
             axis.text.y = element_text(face = "bold", size=12), 
             legend.title = element_text(face = "bold", size=10),
             legend.text = element_text(face = "italic", size=12), 
             legend.key.size = unit(0.5, "cm"),
             panel.grid.minor = element_blank(), 
             panel.background = element_blank())

ggsave("Results/genus_abundance.png", width = unit(15,"cm")) 


# Relative Abundance on Family level

physeq_family_cutoff = preDA(physeq, min.reads = 1000)
glom_or_fam <- tax_glom(physeq_family_cutoff, taxrank ="Family")


gpt_fam <- subset_taxa(glom_or_fam)
gpt_fam <- prune_taxa(names(sort(taxa_sums(glom_or_fam),TRUE)[1:15]), glom_or_fam)
top15_fam <- names(sort(taxa_sums(gpt_fam), decreasing=TRUE))[1:15]
ps.top15_fam <- transform_sample_counts(gpt_fam, function(OTU) 100* OTU/sum(OTU))


vf2=c("peachpuff2","burlywood2","Peachpuff","blanchedalmond", 
      "palegoldenrod","lemonchiffon1", "lightgoldenrod1", "brown1",
      "darkseagreen1",  "seagreen1","darkseagreen3","steelblue3",  
      "steelblue1","brown3","red2","steelblue3",  "steelblue1", "skyblue1",
      "lightsteelblue1", "lightblue2","lightsteelblue", "lightcyan1","coral",
      "sienna1","tan1","orange","brown3","lightsalmon","red2")
col = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
        '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', 
        '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
        '#000075', '#808080', '#ffffff', '#000000', "red", "pink", "cornsilk2")
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
             axis.title.y =element_text(face = "bold", size=11),
             axis.text.x  = element_text(face = "italic", size=12), 
             axis.text.y = element_text(face = "bold", size=12), 
             legend.title = element_text(face = "bold", size=10),
             legend.text = element_text(face = "italic", size=12), 
             legend.key.size = unit(0.5, "cm"),
             panel.grid.minor = element_blank(), 
             panel.background = element_blank())

ggsave("Results/family_abundance_top15.png", width = unit(15,"cm")) 



# I use what is above this line
# _______________________________________________________________________________________



# Do subset on genuslevel! (only the ones we're interested in!)
# We could also subset for species in the same way
physeq_bac_vir <- subset_taxa(physeq, Species == " Porcine circovirus 2" | 
                                Species == "Brachyspira piloscoli" | 
                                Species == " Leptospira interrogans" | 
                                Species == " Bordetella bronchiseptica" | 
                                Species == " Pasteurella multocida" | 
                                Species == " Clostridium difficile" | 
                                Species == " Clostridium tetani" | 
                                Species == " Listeria monocytogenes" | 
                                Species == " Streptococcus suis" | 
                                Species == " Campylobacter jejuni" | 
                                Species == " Glaesserella parasuis" | 
                                Species == " Porcine epidemic diarrhea virus" | 
                                Species == " Clostridium perfringens" | 
                                Species == " Yersinia pseudotuberculosis" | 
                                Species == " Yersinia enterocolitica")
# Missing species: Swine fever, Porcine parvovirus, Pseudorabies virus
sum(sample_sums(physeq_bac_vir))*0.00005 # only use taxa that have at least a 0.00005 fraction of total reads
physeq_bac_vir_cutoff = preDA(physeq_bac_vir, min.reads = 1000) # should we use this number or another one??


physeq_bac_vir_cutoff_df <- psmelt(physeq_bac_vir_cutoff)


ggplot(physeq_bac_vir_cutoff_df, aes(x = Sample, y = Abundance, fill = Species)) + 
  theme_bw() +   
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle=90, size=6))

ggsave("Results/kaijubarplot.png", width = unit(15,"cm"))



# -----------------------------------non-greedy-----------------------------

otutab <- read.csv("non_greedy_kaiju.otu.tab", sep = "\t", row.names = 1, header = TRUE)

otutab_oknames <- otutab %>% 
  rename_at(vars(starts_with("HTA")), funs(str_replace(., "HTAdapter", "Sample")))
otutab_oknames <- otutab_oknames %>% 
  rename_at(vars(starts_with("Sam")), funs(str_replace(., "_FKDL202605956.1a.AK", "_")))
otutab_oknames <- otutab_oknames %>% 
  rename_at(vars(starts_with("Sam")), funs(str_replace(., "_10.*$", "")))
otutab_oknames <- otutab_oknames %>% 
  rename_at(vars(starts_with("Sam")), funs(str_replace(., "_6.*$", "")))
otutab_oknames <- otutab_oknames %>% 
  rename_at(vars(starts_with("Sam")), funs(str_replace(., "_4.*$", "")))
otutab_oknames <- otutab_oknames %>% 
  rename(SampleJapan = DRR025070_1.fq)

# otutab_oknames <- otu_t %>% 
# rename_at(vars(starts_with("DRR")), funs(str_replace(., "DRR025070_1.fq", "SampleJapan")))

colnames(otutab_oknames)

OTU = otu_table(otutab_oknames, taxa_are_rows = TRUE)

taxtab <- read.csv("non_greedy_kaiju.tax.tab", sep = "\t", row.names = 1, header = TRUE)
taxmat = as.matrix(taxtab)
TAX = tax_table(taxmat)

setdiff(rownames (otutab_oknames), rownames(taxmat))

metadata = read.csv("metadata.tsv", sep = "\t", skip=1, header=FALSE)
metadata <- rename(metadata, SampleID=1, Country=2, SampleType=3) #renaming the columns

META = sample_data(metadata)
rownames(META) <-metadata$SampleID
physeq = phyloseq(OTU, TAX, META)
physeq 

saveRDS(physeq, "scrofa.phyloseq_nong.rds")


physeq_df <- psmelt(physeq)
positions=c("Sample8", "Sample9", "Sample10", "Sample11", "Sample12", "Sample13", "Sample14", "Sample15", "Sample16", "Sample17", "Sample18", "Sample19", "Sample20", "Sample21", "SampleJapan", "Sample1", "Sample2", "Sample6", "Sample7", "Sample22", "Sample23")
rawkaijubarplot <- ggplot(physeq_df, aes(x = Sample, y = Abundance, fill = Domain)) + theme_bw() +   geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle=90, size=6)) + scale_x_discrete(limits = positions)
rawkaijubarplot
# dir.create("Results")
ggsave("Results/rawkaijubarplotnong.png") 



# Cutoff species level
physeq_bac_vir_nong <- subset_taxa(physeq, Species == " Porcine circovirus 2" | 
                                     Species == " Leptospira interrogans" | 
                                     Species == " Bordetella bronchiseptica" | 
                                     Species == " Pasteurella multocida" | 
                                     Species == " Clostridium difficile" | 
                                     Species == " Clostridium tetani" | 
                                     Species == " Listeria monocytogenes" | 
                                     Species == " Streptococcus suis" | 
                                     Species == " Campylobacter jejuni" | 
                                     Species == " Glaesserella parasuis" | 
                                     Species == " Porcine epidemic diarrhea virus" | 
                                     Species == " Clostridium perfringens" | 
                                     Species == " Yersinia pseudotuberculosis" | 
                                     Species == " Yersinia enterocolitica")
# Missing species: Swine fever, Porcine parvovirus, Pseudorabies virus


sum(sample_sums(physeq_bac_vir_nong))*0.00005 # only use taxa that ave at least a 0.00005 fraction of total reads
physeq_bac_vir_nong_cutoff = preDA(physeq_bac_vir_nong, min.reads = 1000) # should we use this number or another one??



physeq_bac_vir_nong_cutoff_df <- psmelt(physeq_bac_vir_nong_cutoff)

vf2= c("cadetblue", "cadetblue3", "coral3", "coral", "darkseagreen", 
       "darkseagreen2", "pink3", "pink") 

positions=c("Sample8", "Sample9", "Sample10", "Sample11", "Sample12", 
            "Sample13", "Sample14", "Sample15", "Sample16", "Sample17",
            "Sample18", "Sample19", "Sample20", "Sample21", "SampleJapan", 
            "Sample1", "Sample2", "Sample6", "Sample7", "Sample22", "Sample23")

ggplot(na.omit(physeq_bac_vir_nong_cutoff_df), aes(x = Sample, y = Abundance, fill = Species)) + 
  theme_bw() + geom_bar(stat = "identity") + 
  ylab("Number of reads") + 
  theme(axis.text.x = element_text(angle=90, size=6)) + 
  scale_x_discrete(limits = positions) +
  scale_fill_manual(values = vf2)
  
ggsave("Results/nongkaijubarplot.png", width = unit(15,"cm"))



# Boxplot on nongreedy pathogens grouped by sample type
#-----------------------------------------------------------------------
totalreads <- colSums(otutab_oknames)
pathogens.df <- na.omit(physeq_bac_vir_nong_cutoff_df)

samplelist <- colnames(otutab_oknames)
#adding the total number of reads for each sample. so we can get the precentage abundance
pathogens.df <- pathogens.df %>% mutate(total_reads = case_when(SampleID == "Sample11" ~ totalreads["Sample11"],
                                                                SampleID == "Sample12" ~ totalreads["Sample12"],
                                                                SampleID == "Sample13" ~ totalreads["Sample13"],
                                                                SampleID == "Sample14" ~ totalreads["Sample14"],
                                                                SampleID == "Sample15" ~ totalreads["Sample15"],
                                                                SampleID == "Sample16" ~ totalreads["Sample16"],
                                                                SampleID == "Sample17" ~ totalreads["Sample17"],
                                                                SampleID == "Sample18" ~ totalreads["Sample18"],
                                                                SampleID == "Sample19" ~ totalreads["Sample19"],
                                                                SampleID == "Sample20" ~ totalreads["Sample20"],
                                                                SampleID == "Sample21" ~ totalreads["Sample21"],
                                                                SampleID == "Sample22" ~ totalreads["Sample22"],
                                                                SampleID == "Sample23" ~ totalreads["Sample23"],
                                                                SampleID == "Sample1" ~ totalreads["Sample1"],
                                                                SampleID == "Sample6" ~ totalreads["Sample6"],
                                                                SampleID == "Sample2" ~ totalreads["Sample2"],
                                                                SampleID == "Sample7" ~ totalreads["Sample7"],
                                                                SampleID == "Sample8" ~ totalreads["Sample8"],
                                                                SampleID == "Sample9" ~ totalreads["Sample9"],
                                                                SampleID == "Sample10" ~ totalreads["Sample10"],
                                                                SampleID == "SampleJapan" ~ totalreads["SampleJapan"]
                                                                ))

pathogens.df <- pathogens.df %>% mutate(percent_abundance = Abundance/total_reads) 
names_boxplot = c("C. jejuni", "C. perfringens", "C. tetani", "C. parasuis", "L. monocytogenes", "P. multocida", "S. suis", "Y. enterocolitica")

pathogens.toplot <- pathogens.df %>%  
  dplyr::select(., "Sample", "percent_abundance", "Country", "SampleType", "Species") %>%
  dplyr::filter(!Sample == "Sample1") %>% 
  mutate(SampleType=replace(SampleType, SampleType =="wild boar feces", "Wild life feces")) 
  
ggplot(pathogens.toplot, mapping = aes(x = percent_abundance,  fill = SampleType)) +
  geom_boxplot() + 
  facet_wrap( ~ Species, ncol=4, n=2) +
  scale_fill_manual(values=c("darkseagreen","lightblue3")) +
  coord_flip() +
  xlab("Abundance in %") +
  theme_set(theme_bw()) +
  theme(axis.ticks.x=element_blank(), 
        axis.text.x=element_blank()) 
  


Camp.jejuni <- pathogens.toplot %>% dplyr::filter(Species == " Campylobacter jejuni")
Camp.jejuni.p <- pairwise.wilcox.test(Camp.jejuni$percent_abundance, Camp.jejuni$SampleType, p.adjust.method = "bonferroni") 

clos.perf <- pathogens.toplot %>% dplyr::filter(Species == " Clostridium perfringens")
clos.perf.p <- pairwise.wilcox.test(clos.perf$percent_abundance, clos.perf$SampleType, p.adjust.method = "bonferroni") 

clos.tetani <- pathogens.toplot %>% dplyr::filter(Species == " Clostridium tetani")
clos.tetani.p <- pairwise.wilcox.test(clos.tetani$percent_abundance, clos.tetani$SampleType, p.adjust.method = "bonferroni") 

Glae.para <- pathogens.toplot %>% dplyr::filter(Species == " Glaesserella parasuis")
Glae.para.p <- pairwise.wilcox.test(Glae.para$percent_abundance, Glae.para$SampleType, p.adjust.method = "bonferroni") 

List.mono <- pathogens.toplot %>% dplyr::filter(Species == " Listeria monocytogenes")
List.mono.p <- pairwise.wilcox.test(List.mono$percent_abundance, List.mono$SampleType, p.adjust.method = "bonferroni") 

Past.multo <- pathogens.toplot %>% dplyr::filter(Species == " Pasteurella multocida")
Past.multo.p <- pairwise.wilcox.test(Past.multo$percent_abundance, Past.multo$SampleType, p.adjust.method = "bonferroni") 

Strep.suis <- pathogens.toplot %>% dplyr::filter(Species == " Streptococcus suis")
Strep.suis.p <- pairwise.wilcox.test(Strep.suis$percent_abundance, Strep.suis$SampleType, p.adjust.method = "bonferroni") 

Yers.entero <- pathogens.toplot %>% dplyr::filter(Species == " Yersinia enterocolitica")
Yers.entero.p <- pairwise.wilcox.test(Yers.entero$percent_abundance, Yers.entero$SampleType, p.adjust.method = "bonferroni") 

p.values.path <- rbind(Camp.jejuni.p$p.value, clos.perf.p$p.value, clos.tetani.p$p.value, Glae.para.p$p.value, List.mono.p$p.value,
                       Past.multo.p$p.value, Strep.suis.p$p.value, Yers.entero.p$p.value)
 
rownames(p.values.path) <- c(unique(pathogens.toplot$Species))

p.values.path

### trying to do boxplot on the entire samples:
# glom_domain <- tax_glom(physeq_genus_cutoff, taxrank ="Domain")
# domainlevel.df <- psmelt(glom_domain)
# domainlevel.df <- na.omit(domainlevel.df)
# 
# #takes a while to run.. still waiting
# domainlevel.df %>%  
#   dplyr::select(., "Sample", "Abundance", "Country", "SampleType", "Domain") %>%
#   dplyr::filter(!Sample == "SampleJapan") %>%
#   mutate(SampleType=replace(SampleType, SampleType =="wild boar feces", "Wild life feces")) %>% 
#   ggplot(., mapping = aes(x = Abundance,  fill = SampleType)) +
#   geom_boxplot() + 
#   facet_wrap( ~ Domain) +
#   coord_flip() +
#   xlab("Abundance in reads") +
#   theme_set(theme_bw()) +
#   theme(axis.ticks.x=element_blank(), 
#         axis.text.x=element_blank())
 


# Continue with Beta-Diversity analysis
#-----------------Setup------------------------------------------------------

# Zero correction
physeq_zc <- transform_sample_counts(physeq_bac_vir_nong_cutoff, function(y) sapply(y, function(x) ifelse(x==0, 1, (1-(sum(y==0)*1)/sum(y))*x)))

# log transformation
physeq_clr <- transform_sample_counts(physeq_zc, function(x) log(x/exp(mean(log(x)))))

# PCA -> reduce the dimensionality of the data
ord_clr <- ordinate(physeq_clr, "RDA")
plot_scree(ord_clr) +
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
ggsave("Results/kaijuclrscree.png")

# Plot PC1 & PC2
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)

# In der n?chsten Zeile: was muss ich da tun?
plot_ordination(physeq, ord_clr, type="samples", color="Country") +
  geom_point(size = 6) +
  geom_text(aes(label=Day), colour="black")
ggsave("Results/kaijuclrPCA1x1.png")

# Continue with distance matrix


###________________Differential abundance___________________
# the number of different OTUs with signifigance. 
library(phyloseq)
library(DAtest)
library(DESeq2)
scrofaPhy = readRDS("scrofa.phyloseq_nong.rds")

# phy_genus <- tax_glom(scrofaPhy, "Genus")

treatdds <- phyloseq_to_deseq2(phy_genus, ~ SampleType)
treatdds <- DESeq(treatdds)

res = results(treatdds, alpha = 0.05)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy_genus)[rownames(sigtab), ], "matrix"))

head(sigtab)
dim(sigtab)




