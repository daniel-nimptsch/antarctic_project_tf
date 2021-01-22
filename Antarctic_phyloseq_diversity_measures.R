################################################################################
#
# Antarctic_phyloseq_diversity_measures.R
#
# For the given tables create phyloseq objects to plot different diversity
# indices 
#
# Written by Daniel Nimptsch 
#
################################################################################

#### Packages ####
library(RColorBrewer)
library(ggplot2)
library(phyloseq)
library(readr)
library(tidyr)

#### Functions ####
# Give the otu_mat and the taxa_mat row names and convert the otu_mat to numeric
format_otu_mat = function(otu_mat) {
  rownames(otu_mat) = otu_mat[,1]
  otu_mat = otu_mat[,-1]
  otu_mat = type.convert(otu_mat)
  return(otu_mat)
}

format_tax_mat = function(tax_mat, name) {
  rownames(tax_mat) = tax_mat[,1]
  tax_mat = tax_mat[,-1]
  # tax_mat[,1] = name
  colnames(tax_mat) = c("taxgroup", "taxa")
  return(tax_mat)
}
  
#### Source custom functions ####
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
setwd(paste(own_cloud_dir, "/programms_daniel/Pipeline/r_pipeline_statistics", sep = ""))
source("pipeline_statistics_custom_phyloseq_functions.R")

#### Working Directory ####
setwd(paste(own_cloud_dir, 
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/final_tables_for_grafics/alpha_diversity_indices", 
            sep = ""))
list.files()

final_allOTUs = read.csv(file = list.files(pattern = "final_all1003OTUs"), header = TRUE, sep = "\t")
rownames(final_allOTUs) = final_allOTUs[,1]
final_allOTUs = final_allOTUs[,-1]
final_allOTUs[is.na(final_allOTUs)] = 0

#### Color ####
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#### Selection ####
# Names from the taxgroup tables
taxgroups = list.files(pattern = "allOTUs_with")

# Manual selection
t = 1

# Automatic
for (t in 1:4) {
  name = taxgroups[t]
  # final_table = read.csv(name, sep = "\t", header = TRUE)
  name = strsplit(strsplit(name, split = ".csv")[[1]][1], split = "_")[[1]][3]
  otu_mat = list.files(pattern = paste("allOTUs_withOnlySchF_", name, sep = ""))
  otu_mat = read.csv(otu_mat, sep = "\t", header = TRUE)
  tax_mat = list.files(pattern = paste("tax_mat_", name, sep = ""))
  tax_mat = read.csv(tax_mat, sep = "\t", header = TRUE)
  
    
  # otu_mat = as.matrix(final_table[,c(1,7:12)])
  # tax_mat = as.matrix(final_table[,c(1,3)])
  otu_mat = format_otu_mat(otu_mat)
  tax_mat = format_tax_mat(tax_mat, name)
  otu_mat = cbind(otu_mat, rep(NA, nrow(otu_mat)))
  otu_mat[,] = 0
  colnames(otu_mat) = colnames(final_allOTUs)
  for (i in 1:nrow(otu_mat)) {
    ind = which(rownames(final_allOTUs) == rownames(otu_mat)[i])
    if (length(ind) == 0) {print(rownames(otu_mat)[i])}
    otu_mat[i,] = final_allOTUs[ind,]
  }
  
  otu_mat = otu_mat[,c("AM31", "AM09", "AM914", "AM06", "AM614", "AS14", "AS15", "SchF")]

  # Create a phyloseq object from the otu_mat and the tax_mat
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(tax_mat))
  physeq = phyloseq(OTU, TAX)
  
  #### Diversity Measures ####
  physeq.pruned  = create_physeq_obj_prunded_diversity(physeq)
  richness = estimate_richness(physeq)
  write.table(richness, paste(name, "_diversity_measurements_merged.csv", sep = ""), sep = "\t", row.names = TRUE, col.names = NA)
  alpha_meas = c("Shannon","InvSimpson", "Observed")
  title = paste("Alpha Diversity Plots", "for", name)
  
  #### Final diversity table for summarized plot ####
  div_table = richness[,c(1,6,8)]
  taxa_name_mat = matrix(data = name, ncol = 2, nrow = nrow(richness))
  colnames(taxa_name_mat) = c("Taxa", "Sample")
  taxa_name_mat[,2] = rownames(div_table)
  div_table = cbind(taxa_name_mat, div_table)
  
  #### After all taxgroups ####
  if (t == 1) {
    final_div_table = div_table
  } else {
    final_div_table = rbind(final_div_table, div_table)
  }
}

# Special table for summarized plot for all taxgroups 
value = final_div_table[,c(3,4,5)]
for (i in 1:3) {
  var_final_div_table = final_div_table[,-c(3,4,5)]
  var_value = matrix(data = value[,i], nrow = nrow(final_div_table), ncol = 2)
  var_value[,2] = colnames(value)[i]
  colnames(var_value) = c("Measure", "Index")
  var_final_div_table = cbind(var_final_div_table, var_value)
  if (i == 1) {
    ultra_final_div_table = var_final_div_table
  } else {
    ultra_final_div_table = rbind(ultra_final_div_table, var_final_div_table)
  }
}
ultra_final_div_table = type_convert(ultra_final_div_table)

#### Color ####
# Color Palette
col_vector = brewer.pal(n = 12, name = 'Paired')
display.brewer.pal(n = 12, name = 'Paired')
col_vector[3] = col_vector[4]
col_vector[4] = col_vector[8]

#### Plots ####
# Sample order for the x-axis
sample_order = c("AM31", "AM09", "AM914", "AM06", "AM614", "AS14", "AS15", "SchF")

#### Save ####
#### Diversity Plots ####
pdf(file =  paste(name, "_diversity_plots_merged.pdf", sep = ""), width = 8.2, height = 7.3)
  p <- plot_richness(physeq.pruned, x = "samples", color = "samples", title = title, measures = alpha_meas)
  p$data$samples <- factor(p$data$samples, levels = sample_order)
  p + geom_point(size = 3)
dev.off()

ultra_final_div_table$Index <- factor(ultra_final_div_table$Index, levels = c("Observed", "Shannon", "InvSimpson"))
pdf(file =  "all_taxgroups_diversity_plots_merged.pdf", width = 10.2, height = 6.3)
  p = ggplot(ultra_final_div_table, aes(Sample, Measure))
  p = p + scale_x_discrete(labels = c("AM31-13", "AM09-13", "AM09-14", "AM06-13", "AM06-14", "AS14-14", "AS15-14", "SchF")) 
  p = p + ggtitle("Alpha Diversity Plots for the different taxgroups")
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12))
  p$data$Sample = factor(p$data$Sample, levels = sample_order)
  p = p + facet_wrap(~Index, scales = "free")
  p = p + geom_boxplot(alpha = 0.7, coef = 1.5)
  p = p + geom_point(aes(y = Measure, colour = Taxa), size = 3) + scale_color_manual(values = col_vector)
  p
dev.off()

pdf(file =  "all_taxgroups_diversity_plots_merged_themeBW.pdf", width = 10.2, height = 6.3)
p = ggplot(ultra_final_div_table, aes(Sample, Measure))
p = p + theme_bw()
p = p + scale_x_discrete(labels = c("AM31-13", "AM09-13", "AM09-14", "AM06-13", "AM06-14", "AS14-14", "AS15-14", "SchF")) 
p = p + ggtitle("Alpha Diversity Plots for the different taxgroups")
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12))
p$data$Sample = factor(p$data$Sample, levels = sample_order)
p = p + facet_wrap(~Index, scales = "free")
p = p + geom_boxplot(alpha = 0.7, coef = 1.5)
p = p + geom_point(aes(y = Measure, colour = Taxa), size = 3) + scale_color_manual(values = col_vector)
p
dev.off()