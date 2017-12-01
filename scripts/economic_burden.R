"""
Created in Jan 2017

@author: Alexander Hauser <alexshauser@gmail.com>

Create a stacked box plot to visualise the economic burden
"""

library(dplyr)
library(streamgraph)
library(ggplot2)

font_import(pattern="[A/a]rial")
getPalette = colorRampPalette(brewer.pal(10, "Set1"))

# Economic burden
burden <- read.csv("../data/170829_economic_burden.csv", sep = ",")
head(burden)
burden[, 'variable']  <- as.character(burden[, 'variable'])
burden[, 'value']  <- as.numeric(burden[, 'value'])/1e6

pp <- ggplot(burden, aes(x=variable, y=value)) +
  geom_bar(aes(fill = section), stat="identity") +
  scale_fill_manual(values=getPalette(15)) +
  scale_x_discrete(limits=c('k_homo','k','p_homo','p'),labels=c("homo-\nzygous","all\nvariants","homo-\nzygous","all\nvariants")) +
  labs(x="category",
       y="economic burden on the NHS per year (million GBP)",
       fill = " ") +
  theme_minimal(base_size = 15)
pp
ggsave(file="../figures/170421_economic_burden.eps", plot=pp, width=8, height=6)
