################################################################################
#
# Antarctic_rarefaction_samples.R
#
#
# Written by Daniel Nimptsch 
#
################################################################################

#### Packages ####
library(RColorBrewer)
library(readr)
library(vegan)
library(tidyr)

#### Source custom functions ####
source("pipeline_statistics_custom_phyloseq_functions.R")

#### Working directory ####
setwd("data/rarefaction/")
list.files()

#### Color ####
# Color Palette
col_vector = brewer.pal(n = 8, name = 'Paired')

#### Rarefaction curve ####
final_table = list.files(pattern = "1003")
final_table = read.csv(final_table, sep = "\t", header = TRUE, row.names = 1)
legend = colnames(final_table)
final_table = create_otu_mat_spread(final_table)
raremax = min(rowSums(final_table))
filename = paste("Rarefaction_curve_samples",".pdf", sep = "")

# Plot
pdf(file = filename, width = 11, height = 8)
  par(mar = c(5.1, 5.1, 4.1, 9), xpd = FALSE)
  p = rarecurve(final_table, step = 20, sample = raremax, col = col_vector, lwd = 1.5, ylab = "OTUs", label = FALSE)
  par(xpd = TRUE)
  title("Rarefaction curve of the OTU reads from the samples")
  legend("right", inset = c(-0.16,0), legend = legend, col = col_vector, cex = 0.8, lwd = 2, title = "samples:")
dev.off()

# Alternative
filename = paste("Rarefaction_curve_samples_alternative",".pdf", sep = "")
pdf(file = filename, width = 11, height = 8)
  par(mar = c(5.1, 5.1, 4.1, 9), xpd = FALSE)
  p = rarecurve(final_table, step = 20, col = col_vector, lwd = 1.5, ylab = "OTUs", label = FALSE)
  par(xpd = TRUE)
  title("Rarefaction curve of the OTU reads from the samples")
  legend("right", inset = c(-0.16,0), legend = legend, col = col_vector, cex = 0.8, lwd = 2, title = "samples:")
dev.off()

setwd("../../")