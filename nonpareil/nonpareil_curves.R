rm(list=ls())

library(tidyverse)
library(broom) 
library(Nonpareil)
# setwd("~/Documents/DTU/thesis/thesis_code/nonpareil/")


##plotting the nonpareil curves:
samples <- read.table('nonpareil_table_standard.csv', sep=',', header=TRUE, as.is=TRUE);

attach(samples)
nps <- Nonpareil.set(file, col = col, labels = name, plot.opts = list(plot.observed = FALSE), plot=F)

plot.Nonpareil.Set(x = nps, col = col, labels = name, legend.opts = F, 
                 plot.observed = FALSE  )  
detach(samples)

#generating nonpareil_indices - for violin plot on coverage 
nonpareil_indencies <- print(nps)
nonpareil_indencies <- data.frame(nonpareil_indencies)

ggplot(nonpareil_indencies,aes(x=factor(0),C)) + 
  geom_violin(fill = "white") + 
  geom_boxplot() + 
  ylab("Coverage")+
  theme(axis.title.x=element_blank(), 
        text = element_text(size=50))
ggsave("coverage_violin.png")

nonpareil_indencies %>% summarise(max = max(C), n = n())

