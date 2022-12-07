################################################################################
#
# Antarctic_dotplot_bitscore_summarized_clone.R
#
# Script used to plot the normalized-bitscore values derived from the BLAST query.
# 
#
# Written by Daniel Nimptsch 
#
################################################################################

#### Packages ####
library(RColorBrewer)
library(tidyverse)

#### Working Directory ####
setwd("data/dotplots/")
file = list.files()

# Load the table for the dotplot
table = read.csv(file = file, sep = ";", header = TRUE)
table = table[,c(1:5)]
for (i in 1:nrow(table)) {
  if (table[i,4] == "") {
    table[i,4] = "no_clone"
  } else {
    table[i,4] = "clone"
  }
}
# Detect rows with empty entries and delete them
empty = which(table[,3] == "")
table = table[-empty,]

colnames(table) = c("Taxonomic_name", "Normalized_bitscore", "OTU_ID", "Clone","Taxgroup")
table$Normalized_bitscore = as.numeric(gsub(",", "\\.", table$Normalized_bitscore))
taxonnames = unique(table$Taxonomic_name)

#### Color ####
# Color Palette
col_vector = brewer.pal(n = 12, name = 'Paired')
col_vector[3] = col_vector[4]
col_vector[4] = col_vector[8]
col_vector = col_vector[1:4]

# Save the plot to pdf
# breaks for the y-axis
breaks = seq(0, 1.85, 0.1)
# limits for the y-axis
limits = c(0.1, 1.85)
title = "Dotplot of the normalized bitscores from the different genera"
# Reorder the x-axis labels to the order specified by the .csv
table$Taxonomic_name <- factor(table$Taxonomic_name, levels = unique(table$Taxonomic_name))

# Dotplot
pdf(file = "Summary_dotplot_jitter_with_clones.pdf", width = 12.2, height = 7.3)
  p = ggplot(table, aes(x = Taxonomic_name, y = Normalized_bitscore))
  p = p + facet_grid(cols = vars(Taxgroup), space = "free_x", scales = "free_x")
  p = p + theme(strip.text.x = element_text(size = 6))
  p = p + scale_y_continuous(limits = limits, breaks = breaks, minor_breaks = FALSE, 
                             sec.axis = sec_axis(~ ., breaks = breaks))
  p = p + geom_jitter(data = table[table$Clone == "no_clone",], width = 0.1, shape = 21, colour = "black", 
                      stroke = 0.2, aes(fill = Taxgroup)) + guides(fill = FALSE)
  p = p + geom_jitter(data = table[table$Clone == "clone",], width = 0.1, shape = 23, colour = "black", 
                      stroke = 0.2, fill = "yellow") + guides(fill = FALSE)
  p = p + scale_fill_manual(values = col_vector)
  p = p + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7))
  p = p + theme(panel.grid.major = element_line(size = 0.2))
  p = p + theme(axis.title.y = element_text(vjust = 2.5))
  p = p + ggtitle(title)
  p
dev.off()

pdf(file = "Summary_dotplot_jitter_with_clones_themeBW.pdf", width = 12.2, height = 7.3)
  p = ggplot(table, aes(x = Taxonomic_name, y = Normalized_bitscore))
  p = p + theme_bw()
  p = p + facet_grid(cols = vars(Taxgroup), space = "free_x", scales = "free_x")
  p = p + theme(strip.text.x = element_text(size = 6))
  p = p + scale_y_continuous(limits = limits, breaks = breaks, minor_breaks = FALSE, 
                             sec.axis = sec_axis(~ ., breaks = breaks))
  p = p + geom_jitter(data = table[table$Clone == "no_clone",], width = 0.1, shape = 21, colour = "black", 
                      stroke = 0.2, aes(fill = Taxgroup)) + guides(fill = FALSE)
  p = p + geom_jitter(data = table[table$Clone == "clone",], width = 0.1, shape = 23, colour = "black", 
                      stroke = 0.2, fill = "yellow") + guides(fill = FALSE)
  p = p + scale_fill_manual(values = col_vector)
  p = p + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7))
  p = p + theme(panel.grid.major = element_line(size = 0.2))
  p = p + theme(axis.title.y = element_text(vjust = 2.5))
  p = p + ggtitle(title)
  p
dev.off()

# Dotplot + Boxplot
pdf(file = "Summary_dotplot_and_boxplot_jitter_with_clones.pdf", width = 12.2, height = 7.3)
  p = ggplot(table, aes(x = Taxonomic_name, y = Normalized_bitscore))
  p = p + facet_grid(cols = vars(Taxgroup), space = "free_x", scales = "free_x")
  p = p + theme(strip.text.x = element_text(size = 6))
  p = p + scale_y_continuous(limits = limits, breaks = breaks, minor_breaks = FALSE, 
                             sec.axis = sec_axis(~ ., breaks = breaks))
  p = p + geom_boxplot(alpha = 0.4, lwd = 0.2, outlier.alpha = 0)
  p = p + geom_jitter(data = table[table$Clone == "no_clone",], width = 0.1, shape = 21, colour = "black", 
                      stroke = 0.2, aes(fill = Taxgroup)) + guides(fill = FALSE)
  p = p + geom_jitter(data = table[table$Clone == "clone",], width = 0.1, shape = 23, colour = "black", 
                      stroke = 0.2, fill = "yellow") + guides(fill = FALSE)
  p = p + scale_fill_manual(values = col_vector)
  p = p + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7))
  p = p + theme(panel.grid.major = element_line(size = 0.2))
  p = p + theme(axis.title.y = element_text(vjust = 2.5))
  p = p + ggtitle(title)
  p
dev.off()

pdf(file = "Summary_dotplot_and_boxplot_jitter_with_clones_themeBW.pdf", width = 12.2, height = 7.3)
  p = ggplot(table, aes(x = Taxonomic_name, y = Normalized_bitscore))
  p = p + theme_bw()
  p = p + facet_grid(cols = vars(Taxgroup), space = "free_x", scales = "free_x")
  p = p + theme(strip.text.x = element_text(size = 6))
  p = p + scale_y_continuous(limits = limits, breaks = breaks, minor_breaks = FALSE, 
                             sec.axis = sec_axis(~ ., breaks = breaks))
  p = p + geom_boxplot(alpha = 0.7, lwd = 0.2, outlier.alpha = 0)
  p = p + geom_jitter(data = table[table$Clone == "no_clone",], width = 0.1, shape = 21, colour = "black", 
                      stroke = 0.2, aes(fill = Taxgroup)) + guides(fill = FALSE)
  p = p + geom_jitter(data = table[table$Clone == "clone",], width = 0.1, shape = 23, colour = "black", 
                      stroke = 0.2, fill = "yellow") + guides(fill = FALSE)
  p = p + scale_fill_manual(values = col_vector)
  p = p + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7))
  p = p + theme(panel.grid.major = element_line(size = 0.2))
  p = p + theme(axis.title.y = element_text(vjust = 2.5))
  p = p + ggtitle(title)
  p
dev.off()

setwd("../../")