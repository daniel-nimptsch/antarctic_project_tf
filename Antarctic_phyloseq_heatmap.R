################################################################################
#
# Antarctic_phyloseq_heatmap.R
#
# For the given tables create phyloseq objects to plot heatmaps of the top 20
# OTUs from the Antarctic project
#
# Written by Daniel Nimptsch 
#
################################################################################

#### Packages ####
library(phyloseq)
library(tidyverse)

#### Functions ####
# Create a phyloseq object with a given taxonomy table and a OTU/ASV-table (as matrix)
create_physeq_obj = function(otu_mat, tax_mat) {
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  physeq = phyloseq(OTU, TAX)
  return(physeq)
}

# Edit of the plot_heatmap function
source("Antarctic_phyloseq_heatmap_custom_functions.R")

#### Working directory ####
setwd("data/heatmaps/")

################################################################################

#### Selection ####
for_heatmaps = list.files(pattern = "forheatmap.csv")

# Which one?
i = 4
name = for_heatmaps[i]

# Read the table and generate a otu_mat and a tax_mat
final_table = read.csv(name, sep = ";", header = TRUE)
final_table[,2] = gsub(",", ".", final_table[,2])

# Remove ".csv" from the name
name = strsplit(name, split = ".csv")[[1]][1]

otu_mat = as.matrix(final_table[,c(1,5:10)])
tax_mat = as.matrix(final_table[,c(1,3)])

colnames(tax_mat) = c("OTU_ID", "species_taxonomy")
rownames(tax_mat) = tax_mat[,1]
tax_mat_right = tax_mat

# Change the taxa string to a more complex: OTU_ID_normalized_bitscore_taxonmy
for (y in 1:nrow(final_table)) {
  tax_mat[y,2] = paste(final_table[y,1], final_table[y,2], final_table[y,3], sep = "_")
}

# Right alignment version
for (y in 1:nrow(final_table)) {
  tax_mat_right[y,2] = paste(final_table[y,3], final_table[y,2], final_table[y,1], sep = "_")
}

# Give the otu_mat and the taxa_mat row names and convert the otu_mat to numeric
rownames(otu_mat) = otu_mat[,1]

otu_mat = otu_mat[,-1]
otu_mat = type.convert(otu_mat)

# Create a phyloseq object from the otu_mat and the tax_mat
physeq = create_physeq_obj(otu_mat, tax_mat)
physeq_right = create_physeq_obj(otu_mat, tax_mat_right)
sample_oder = c("AM31", "AM09", "AM06", "AS14", "AS15", "SchF")
taxa_order = rev(as.vector(tax_mat[,1]))
taxname = strsplit(name, split = "_")[[1]][1]
title = paste("Heatmap of the OTU-abundacies from the top 20 ", taxname, " Taxas", sep = "")

################################################################################
# Print and save the heatmap to pdf
# Custom version with edited source code -> y-axis to the right
pdf(file = paste(name, "_dummy_right.pdf", sep = ""), width = 5.6, height = 8.3)
  p = plot_heatmap_ypos_right(physeq_right, low = "#000033", high = "#CCFF66", 
                              taxa.label = "species_taxonomy", taxa.order = taxa_order,
                              title = title, sample.order = sample_oder)
  p = p + theme(axis.text.x = element_text(size = 10, angle = -45), 
                axis.text.y = element_text(size = 10, angle = 0),
                legend.position = "bottom",
                legend.text = element_text(size = 7),
                plot.title = element_text(size = 11, hjust = 1))
  p
dev.off()

# Default version -> y-axis to the left
pdf(file = paste(name, "_dummy_left.pdf", sep = ""), width = 5.6, height = 8.3)
  p = plot_heatmap(physeq, low = "#000033", high = "#CCFF66", 
                   taxa.label = "species_taxonomy", taxa.order = taxa_order,
                   title = title, sample.order = sample_oder)
  p = p + theme(axis.text.x = element_text(size = 10, angle = -45), 
                axis.text.y = element_text(size = 10, angle = 0),
                legend.position = "bottom",
                legend.text = element_text(size = 7),
                plot.title = element_text(size = 11, hjust = 1))
  p
dev.off()

setwd("../../")
