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
  tax_mat[,1] = name
  colnames(tax_mat) = c("taxgroup", "taxa")
  return(tax_mat)
}
  
#### Source custom functions ####
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
setwd(paste(own_cloud_dir, "/programms_daniel/Pipeline/r_pipeline_statistics", sep = ""))
source("pipeline_statistics_custom_phyloseq_functions.R")

#### Working Directory ####
setwd(paste(own_cloud_dir, "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/final_tables_for_grafics/alpha_diversity_indices", sep = ""))
list.files()

#### Color ####
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#### Selection ####
# Names from the taxgroup tables
taxgroups = list.files(pattern = "final_allOTUs")

# Manual selection
# t = 1

# Automatic
for (t in 1:4) {
  name = taxgroups[t]
  final_table = read.csv(name, sep = "\t", header = TRUE)
  name = strsplit(strsplit(name, split = ".csv")[[1]][1], split = "_")[[1]][3]
  otu_mat = as.matrix(final_table[,c(1,7:12)])
  tax_mat = as.matrix(final_table[,c(1,3)])
  otu_mat = format_otu_mat(otu_mat)
  tax_mat = format_tax_mat(tax_mat, name)
  otu_mat[is.na(otu_mat)] = 0

  #### Add the only_SchF OTUs ####
  only_SchF_table = read.csv("/home/pipeline/bigdata/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/only_SchF/OtuMatrix_onlySchF.csv", sep = "\t", header = TRUE)
  only_SchF_table = only_SchF_table[,-c(2:15)]
  only_SchF_table[,2] = only_SchF_table[,2] + only_SchF_table[,3]
  only_SchF_table = only_SchF_table[,-c(3,4)]
  colnames(only_SchF_table) = c("OTU_ID", "only_SchF", "taxgroup")
  # merge the tables
  only_SchF_table = only_SchF_table[which(only_SchF_table$taxgroup == name),]
  merge_only_SchF = matrix(data = 0, nrow = nrow(otu_mat) + nrow(only_SchF_table), ncol = ncol(otu_mat) + 1)
  colnames(merge_only_SchF) = append(colnames(otu_mat), "only_SchF")
  row.names(merge_only_SchF) = append(rownames(otu_mat), only_SchF_table[,1])
  for (i in 1:nrow(merge_only_SchF)) {
    if (length(which(row.names(merge_only_SchF)[i] == row.names(otu_mat))) != 0) {
      ind = which(row.names(merge_only_SchF)[i] == row.names(otu_mat))
      merge_only_SchF[i,c(1:6)] = otu_mat[ind,]
    } else if (length(which(row.names(merge_only_SchF)[i] == only_SchF_table[,1])) != 0) {
      ind = which(row.names(merge_only_SchF)[i] == only_SchF_table[,1])
      merge_only_SchF[i,7] = only_SchF_table[ind,2]
    }
  }
  merge_only_SchF[,7] = merge_only_SchF[,6] + merge_only_SchF[,7]
  otu_mat = merge_only_SchF

  # merge the tax_mat
  only_SchF_table_tax = only_SchF_table[,c(2,3)]
  rownames(only_SchF_table_tax) = only_SchF_table[,1]
  only_SchF_table_tax[,1] = only_SchF_table_tax[,2]
  colnames(only_SchF_table_tax) = c("taxgroup", "taxa")
  only_SchF_table_tax = as.matrix(only_SchF_table_tax)
  tax_mat = rbind(tax_mat, only_SchF_table_tax)
  
  # Create a phyloseq object from the otu_mat and the tax_mat
  physeq = create_physeq_obj(otu_mat, tax_mat)
  
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
sample_order = c("AM31", "AM09", "AM06", "AS14", "AS15", "SchF", "only_SchF")

#### Save ####
#### Diversity Plots ####
pdf(file =  paste(name, "_diversity_plots_merged.pdf", sep = ""), width = 8.2, height = 7.3)
  p <- plot_richness(physeq.pruned, x = "samples", color = "samples", title = title, measures = alpha_meas)
  p$data$samples <- factor(p$data$samples, levels = sample_order)
  p + geom_point(size = 3)
dev.off()

ultra_final_div_table$Index <- factor(ultra_final_div_table$Index, levels = c("Observed", "Shannon", "InvSimpson"))
pdf(file =  "all_taxgroups_diversity_plots_merged.pdf", width = 10.2, height = 7.3)
  p = ggplot(ultra_final_div_table, aes(Sample, Measure))
  p = p + scale_x_discrete(labels = c("AM31", "AM09", "AM06", "AS14", "AS15", "shared_SchF", "SchF")) 
  p = p + ggtitle("Alpha Diversity Plots for the different taxgroups")
  p = p + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  p$data$Sample = factor(p$data$Sample, levels = sample_order)
  p = p + facet_wrap(~Index, scales = "free")
  p = p + geom_boxplot(alpha = 0, coef = 1.5)
  p = p + geom_point(aes(y = Measure, colour = Taxa), size = 3) + scale_color_manual(values = col_vector)
  p
dev.off()