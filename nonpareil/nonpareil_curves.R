rm(list=ls())

library(tidyverse)
library(broom) 
library(Nonpareil)
setwd("~/Documents/DTU/thesis/thesis_code/nonpareil/")


##plotting the nonpareil curves:
samples <- read.table('nonpareil_table_standard.csv', sep=',', header=TRUE, as.is=TRUE);
# samples <- read.table('nonpareil_table_60.csv', sep=',', header=TRUE, as.is=TRUE);

attach(samples)
nps <- Nonpareil.set(file, col=col, labels=name, plot.opts=list(plot.observed=FALSE), plot=F)

pplot.Nonpareil.Set(x = nps, col = col, labels = name, legend.opts = F, 
                   plot.observed = FALSE) 
# Nonpareil.legend(nps, x="bottomleft", cex=0.45)
detach(samples)

#generating nonpareil_indices - for violin plot on coverege 
nonpareil_indencies <- print(nps)
nonpareil_indencies <- data.frame(nonpareil_indencies)

ggplot(nonpareil_indencies,aes(x=factor(0),C)) + 
  geom_violin(fill = "white") + 
  geom_boxplot() + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank())
ggsave("coverage_violin.png")
