# nonpareil for project 
rm(list=ls())

library(tidyverse)
library(broom) 
library(Nonpareil)
setwd("~/Documents/DTU/Semester5/metagenomicsMicrobiomeAnalysis/project/")
metadata <- read_tsv("metadata.tsv", skip=1, col_names=FALSE)
dataset <- read_csv("nonpareil_indices.csv", skip=1, col_names=FALSE) 
datasetJapan <- read_csv("nonpareil_indices_japan.csv", skip=1, col_names=FALSE) 

head(dataset)
datasetJapan <- rename(datasetJapan, Sample = 1, kappa=2, C=3, LR=4, modelR=5, LRstar=6, diversity=7) #renaming the columns
dataset <- rename(dataset, Sample = 1, kappa=2, C=3, LR=4, modelR=5,LRstar=6,diversity=7) #renaming the columns
dataset <- dataset %>% mutate_at('Sample', funs(str_replace(., "HTAdapter", "Sample")))  
dataset <- dataset %>% mutate_at('Sample', funs(str_replace(., "_FKDL202605956.1a.AK", "_")))
dataset <- dataset %>% mutate_at('Sample', funs(str_replace(., "_10.*$", "")))
dataset <- dataset %>% mutate_at('Sample', funs(str_replace(., "_6.*$", "")))
dataset <- dataset %>% mutate_at('Sample', funs(str_replace(., "_4.*$", "")))
datasetJapan <- datasetJapan %>% mutate_at('Sample', str_replace, "DRR025070",  "SampleJapan")  
data <- rbind(datasetJapan, dataset)

metadata <- rename(metadata, SampleID=1, Country=2, SampleType=3) #renaming the columns
dataset <- full_join(data, metadata, by = c("Sample" = "SampleID")) 

ggplot(dataset,aes(x=factor(0),C)) + 
  geom_violin(fill = "Grey") + 
  geom_boxplot() + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())
# ggsave("coverage.png") 
ggplot(dataset, aes(diversity, C)) + geom_point()
# ggsave("coverage_diversity.png")

#checking the correlation between coverage and diversity. 
lm(C ~ diversity, data = dataset) %>% tidy 
ggplot(dataset, aes(x = diversity, y = C)) + geom_point() + geom_smooth(method = "lm") + annotate("text", x = 19, y = 1, label = "-0.136±0.009\np-val 6.06e-15" ) 

#plotting the diversity against sample type. 
ggplot(data = dataset) + 
  geom_violin(mapping = aes(x =Country , y = diversity, fill = SampleType), trim = FALSE)
# ggsave("sampletypeVSdivViolin.png") 

data_2 <- dataset %>% 
  mutate(SampleType=replace(SampleType, SampleType =="wild boar feces", "Wild life feces"))

pairwise.wilcox.test(data_2$diversity, data_2$SampleType, p.adjust.method = "bonferroni") 



data_2 %>% select("SampleType", "C") %>%  summarise(m =mean(C))
# mean coverage for the two sample types:
coverage_mean <- data_2 %>% select("SampleType", "C") %>% group_by(SampleType) %>% summarise(m =mean(C))

##plotting the nonpareil curves:
samples <- read.table('nonpareil.table.txt', sep='\t', header=TRUE, as.is=TRUE);

attach(samples)
nps <- Nonpareil.set(file, col=col, labels=name, plot.opts=list(plot.observed=FALSE), plot=F)

plot.Nonpareil.Set(x = nps, col = col, labels = name, legend.opts = F, 
                   plot.observed = FALSE)
Nonpareil.legend(nps, x="bottomleft", cex=0.55)
detach(samples)

# #generating nonpareil_indices - for violin plot on coverege based on nps
# nonpareil_indencies <- print(nps)
