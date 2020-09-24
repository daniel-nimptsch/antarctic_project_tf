#### Packages ####
library(RColorBrewer)
library(readr)
library(vegan)
library(tidyr)

#### Source custom functions ####
setwd("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline/r_pipeline_statistics")
source("pipeline_statistics_custom_phyloseq_functions.R")

#### Working directory ####
setwd("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/final_OTUs_NGS_Ant1/processed")
list.files()

#### Load the final_tables ####

# Final table
final_tables = list.files(pattern = ".csv")
row_final_tables = 0
colnames_final_tables = c()
colnames_final_tables[1] = "OTU_ID"
for (i in 1:length(final_tables)) {
  row_final_tables = row_final_tables + nrow(read.csv(final_tables[i], sep = ",", header = TRUE, row.names = 1))
  colnames_final_tables[i+1] = strsplit(strsplit(final_tables[i], split = "_")[[1]][4], split = "\u002Ecsv")[[1]][1]
}
final_otu_table = matrix(data = NA, nrow = row_final_tables, ncol = 5)
colnames(final_otu_table) = colnames_final_tables

for (i in 1:length(final_tables)) {
  final_table = read.csv(final_tables[i], sep = ",", header = TRUE)
  name = strsplit(strsplit(final_tables[i], split = "_")[[1]][4], split = "\u002Ecsv")[[1]][1]
  # generate a otu_mat
  otu_mat = final_table[,c(1,6,7,8,9,10)]
  read_sums = matrix(data = NA, nrow = nrow(otu_mat), ncol = 1)
  # calculate the rowsums
  for (y in 1:nrow(otu_mat)) {
    read_sums[y,1] = sum(otu_mat[y,2:6])
  }
  otu_mat = cbind(otu_mat, read_sums)
  otu_mat = otu_mat[,c(1,7)]
  
  col_ind = which(colnames(final_otu_table) == name)
  for (y in 1:nrow(otu_mat)) {
    ind = which(is.na(final_otu_table[,1]))[1]
    final_otu_table[ind,1] = otu_mat[y,1]
    final_otu_table[ind,col_ind] = otu_mat[y,2]
  }  
}
final_otu_table[is.na(final_otu_table)] = 0

# bak
final_otu_table_bak = final_otu_table
final_otu_table = final_otu_table_bak

#### Color ####
# Color Palette
col_vector = brewer.pal(n = 12, name = 'Paired')
col_vector[3] = col_vector[4]
col_vector[4] = col_vector[8]

# cloro
final_otu_table = as.matrix(final_otu_table[,2])
final_otu_table = create_otu_mat_spread(final_otu_table)
col_vector = col_vector[1]
final_otu_table = matrix(as.numeric(unlist(final_otu_table)),nrow=nrow(final_otu_table))
raremax = min(rowSums(final_otu_table))
legend = colnames_final_tables[2]

# clorotreboux
final_otu_table = final_otu_table[,c(2,3)]
final_otu_table = create_otu_mat_spread(final_otu_table)
col_vector = col_vector[1:2]
final_otu_table = matrix(as.numeric(unlist(final_otu_table)),nrow=nrow(final_otu_table))
raremax = min(rowSums(final_otu_table ))
legend = colnames_final_tables[2:3]

# treboux
final_otu_table = as.matrix(final_otu_table[,3])
final_otu_table = create_otu_mat_spread(final_otu_table)
col_vector = col_vector[2]
final_otu_table = matrix(as.numeric(unlist(final_otu_table)),nrow=nrow(final_otu_table))
raremax = min(rowSums(final_otu_table))
legend = colnames_final_tables[3]

# ulvo
final_otu_table = as.matrix(final_otu_table[,4])
final_otu_table = create_otu_mat_spread(final_otu_table)
col_vector = col_vector[3]
final_otu_table = matrix(as.numeric(unlist(final_otu_table)),nrow=nrow(final_otu_table))
raremax = min(rowSums(final_otu_table))
legend = colnames_final_tables[4]

# xantho
final_otu_table = as.matrix(final_otu_table[,5])
final_otu_table = create_otu_mat_spread(final_otu_table)
col_vector = col_vector[4]
final_otu_table = matrix(as.numeric(unlist(final_otu_table)),nrow=nrow(final_otu_table))
raremax = min(rowSums(final_otu_table))
legend = colnames_final_tables[5]


#########################################################################################
#### Plots ####
if(!dir.exists("R_Statistik")) { dir.create("R_Statistik")}

#### Rarefaction curve ####
pdf(file = "R_Statistik/Rarefaction_curve_taxgroups_ChloroTreboux.pdf", width = 11, height = 8)
par(mar=c(5.1, 5.1, 4.1, 9), xpd = FALSE)
p = rarecurve(final_otu_table, step = 20, sample = raremax, lwd = 2.5, 
              col = col_vector, ylab="OTUs", label=F) 
title("Rarefaction curve of the OTU counts from the taxgroups")
par(xpd=TRUE)
legend("right", inset = c(-0.21,0), legend = legend, cex = 0.6, lwd = 2.5, 
           col = col_vector, title = "Taxgroups:")
dev.off()