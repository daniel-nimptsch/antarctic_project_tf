################################################################################
#
# Antarctic_rarefaction.R
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
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
setwd(paste(own_cloud_dir, "/programms_daniel/Pipeline/r_pipeline_statistics", sep = ""))
source("pipeline_statistics_custom_phyloseq_functions.R")

#### Working directory ####
setwd(paste(own_cloud_dir, "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/final_tables_for_grafics/rarefaction", sep = ""))
list.files()

#### Functions ####
# Function to prepare the data for the rarefaction plot of the corresponing taxgroup
# and to actually draw and save the plot to pdf
plot_rarefaction = function(final_otu_table, col_vector, column_nr) {
  final_otu_table = as.matrix(final_otu_table[,column_nr])
  final_otu_table = create_otu_mat_spread(final_otu_table)
  final_otu_table = matrix(as.numeric(unlist(final_otu_table)),nrow=nrow(final_otu_table))
  raremax = min(rowSums(final_otu_table))
  legend = colnames_final_tables[column_nr]
  #### Rarefaction curve ####
  filename = paste("Rarefaction_curve_taxgroups_", legend,".pdf", sep = "")
  pdf(file = filename, width = 11, height = 8)
  par(mar=c(5.1, 5.1, 4.1, 9), xpd = FALSE)
  p = rarecurve(final_otu_table, step = 20, sample = raremax, lwd = 2.5, 
                col = col_vector, ylab="OTUs", label=F) 
  title("Rarefaction curve of the OTU counts from the taxgroups")
  par(xpd=TRUE)
  legend("right", inset = c(-0.21,0), legend = legend, cex = 0.6, lwd = 2.5, 
         col = col_vector, title = "Taxgroups:")
  dev.off()
}

#### Load the final_tables ####
# Final table
final_tables = list.files(pattern = ".csv")
# Empty variables
row_final_tables = 0
colnames_final_tables = c()
colnames_final_tables[1] = "OTU_ID"
# For all the taxgroups count the total number of OTUs
for (i in 1:length(final_tables)) {
  row_final_tables = row_final_tables + nrow(read.csv(final_tables[i], sep = "\t", header = TRUE, row.names = 1))
  colnames_final_tables[i+1] = strsplit(strsplit(final_tables[i], split = "_")[[1]][3], split = "\\.")[[1]][1]
}
final_otu_table = matrix(data = NA, nrow = row_final_tables, ncol = 5)
colnames(final_otu_table) = colnames_final_tables
# Determine for all the taxgroups the sum of the reads from the OTUs
for (i in 1:length(final_tables)) {
  final_table = read.csv(final_tables[i], sep = "\t", header = TRUE)
  name = strsplit(strsplit(final_tables[i], split = "_")[[1]][3], split = "\\.")[[1]][1]
  # Generate a otu_mat
  otu_mat = final_table[,c(1,7:11)]
  otu_mat[is.na(otu_mat)] = 0
  read_sums = matrix(data = NA, nrow = nrow(otu_mat), ncol = 1)
  # Calculate the rowsums
  for (y in 1:nrow(otu_mat)) {
    read_sums[y,1] = sum(otu_mat[y,2:6])
  }
  otu_mat = cbind(otu_mat, read_sums)
  # otu_mat now only with the sums
  otu_mat = otu_mat[,c(1,7)]
  # Add the row sums to the final table
  col_ind = which(colnames(final_otu_table) == name)
  for (y in 1:nrow(otu_mat)) {
    ind = which(is.na(final_otu_table[,1]))[1]
    final_otu_table[ind,1] = otu_mat[y,1]
    final_otu_table[ind,col_ind] = otu_mat[y,2]
  }  
}

#### Color ####
# Color Palette
col_vector = brewer.pal(n = 12, name = 'Paired')
col_vector[3] = col_vector[4]
col_vector[4] = col_vector[8]

# cloro
plot_rarefaction(final_otu_table, col_vector[1], 2)

# treboux
plot_rarefaction(final_otu_table, col_vector[2], 3)

# ulvo
plot_rarefaction(final_otu_table, col_vector[3], 4)

# xantho
plot_rarefaction(final_otu_table, col_vector[4], 5)
