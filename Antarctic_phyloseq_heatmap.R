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
library(ggplot2)
library(phyloseq)
library(readr)
library(tidyr)

#### Source custom functions ####
setwd("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline/r_pipeline_statistics")
source("pipeline_statistics_custom_phyloseq_functions.R")

#### Working directory ####
setwd("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/heatmaps")

#### Selection ####
for_heatmaps = list.files(pattern = "heatmap_top20.csv")
# special
name = list.files(pattern = "special_heatmap.csv")
# which one?
# i = 4
# name = for_heatmaps[i]

# read the table and generate a otu_mat and a tax_mat
# special: "\t", other ";"
final_table = read.csv(name, sep = "\t", header = TRUE)
final_table[,2] = gsub(",", ".", final_table[,2])
# remove ".csv" from the name
name = strsplit(name, split = ".csv")[[1]][1]
# Special 9, other 10
otu_mat = as.matrix(final_table[,c(1,5:10)])
tax_mat = as.matrix(final_table[,c(1,3)])
# change the taxa string to a more complex: OTU_ID_normalized_bitscore_taxonmy
for (y in 1:nrow(final_table)) {
  tax_mat[y,2] = paste(final_table[y,1], final_table[y,2], final_table[y,3], sep = "_")
}
# Give the otu_mat and the taxa_mat row names and convert the otu_mat to numeric
rownames(otu_mat) = otu_mat[,1]
otu_mat = otu_mat[,-1]
otu_mat = type.convert(otu_mat)
rownames(tax_mat) = tax_mat[,1]
# create a phyloseq object from the otu_mat and the tax_mat
physeq = create_physeq_obj(otu_mat, tax_mat)
sample_oder = c("AM31", "AM09", "AM06", "AS14", "AS15", "SchF")
taxa_order = rev(as.vector(tax_mat[,1]))
# taxa_order = rev(taxa_names(physeq))
# print and save the heatmap to pdf
pdf(file = paste(name, ".pdf", sep = ""), width = 8.2, height = 8.3)
p = plot_heatmap(physeq, low="#000033", high="#CCFF66", taxa.label = "species.taxonomy", taxa.order = taxa_order,
                 title = "Heatmap of the OTU-abundacies", sample.order = sample_oder)
p + theme(axis.text.x = element_text(angle = -45), axis.text.y = element_text(size = 9))
dev.off()