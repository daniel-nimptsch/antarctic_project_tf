################################################################################
# 
# Antarctic_NMDS.R
#
# This script is used to generate NMDS ordination plots from the taxgroups: 
# "Chlorophyceae", "Trebouxiophyceae", "Ulvophyceae" and "Xanthophyceae" with the 
# Hellinger transformation. It saves the plots to pdf and aditionally the OTU- and
# tax-table to csv.
#
# Written by Daniel Nimptsch
#
################################################################################

#### Packages ####
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(stringr)
library(plyr)
library(vegan)

#### Functions ####
# Create a phyloseq object with a given taxonomy table and a OTU/ASV-table (as matrix)
create_physeq_obj = function(otu_mat, tax_mat) {
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  physeq = phyloseq(OTU, TAX)
  return(physeq)
}

#### Working Machine ####
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")

#### Working directory ####
setwd(paste(own_cloud_dir, 
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/final_tables_for_grafics/NDMS"
            , sep = ""))
path = getwd()
list.files()

# Load the table
table = read.csv(list.files(pattern = "fNMDS"), header = TRUE, sep = "\t")
table[is.na(table)] = 0

#### otu-Matrix ####
otu_mat = matrix(data = NA, nrow = nrow(table), ncol = 7)
colnames(otu_mat) = c("OTU_ID", "AS14", "AS15", "AM31", "AM06", "AM09", "SchF")
otu_mat[,1] = table[,1]
for (y in 1:nrow(table)) {
  # AS14
  otu_mat[y,2] = sum(table[y,grep("X14", colnames(table))])
  # AS15
  otu_mat[y,3] = sum(table[y,grep("X15", colnames(table))])
  # AS31
  otu_mat[y,4] = sum(table[y,grep("X31", colnames(table))])
  # AS06  
  otu_mat[y,5] = sum(table[y,grep("X6", colnames(table))])
  # AS09
  otu_mat[y,6] = sum(table[y,grep("X9", colnames(table))])
  # SchF
  otu_mat[y,7] = sum(table[y,grep("SchF", colnames(table))])
}
rownames(otu_mat) = otu_mat[,1]
otu_mat = otu_mat[,-1]
otu_mat = type.convert(otu_mat)

# Exclude SchF
otu_mat = otu_mat[,-6]

#### tax-Matrix ####
tax_mat = matrix(data = NA, nrow = nrow(table), ncol = 2)
tax_mat[,1] = rownames(otu_mat)
rownames(tax_mat) = tax_mat[,1]
tax_mat[,c(1,2)] = table$taxgroup
colnames(tax_mat) = c("taxgroup_1", "taxgroup")
tax_mat[,2] = "other"

#### taxgroups ####
# Use the information from the taxgroup files to complement the tax-table
taxgroup_files = list.files(pattern = "final_allOTUs")
for (i in 1:length(taxgroup_files)) {
  taxgroup = read.csv(taxgroup_files[i], header = TRUE, sep = "\t")
  name = str_remove(taxgroup_files[i], "final_allOTUs_")
  name = str_remove(name, ".csv")
  tax_mat[taxgroup[,1],2] = name
}

# Generate phyloseq objects for the different taxgroups.
# Then apply the hellinger transformation to the OTU-table within the 
# phyloseq object and then do the nmds ordination for the taxgroups individually.
# Finally save the ordination plots to a list (p-list, plot-list)
plist = list()
ord_list = list()
for (i in 1:length(taxgroup_files)) {
  name = str_remove(taxgroup_files[i], "final_allOTUs_")
  name = str_remove(name, ".csv")
  otu_mat_taxgr = otu_mat[which(tax_mat[,2] == name),]
  tax_mat_taxgr = tax_mat[which(tax_mat[,2] == name),]
  
  #### Phyloseq ####
  physeq = create_physeq_obj(otu_mat_taxgr, tax_mat_taxgr)
  Location = as.data.frame(sample_names(physeq))
  colnames(Location) = "Location"
  Location = sample_data(data.frame(Location, row.names = sample_names(physeq)))
  physeq = merge_phyloseq(physeq, Location)
  sample_data(physeq)$taxgroup <- name
  
  #### Hellinger Transformation ###
  physeq.hell = create_physeq_obj(otu_mat_taxgr, tax_mat_taxgr)
  sample_data(physeq.hell) = sample_data(physeq)
  physeq.hell <- transform_sample_counts(physeq, function(x) sqrt(x / sum(x)))
  physeq_ord = ordinate(physeq.hell, "NMDS", "bray")
  Sys.sleep(1)
  plist[[i]] = plot_ordination(physeq.hell, physeq_ord, type = "sample", color = "taxgroup", shape = "Location")
  ord_list[[i]] = physeq_ord
}

# Generate a dataframe from the generated p-list (plot-list)
dataframe = data.frame()
for (i in 1:length(plist)) {
  dataframe = rbind(dataframe, as.data.frame(plist[[i]]$data))
}
row.names(dataframe) = seq(nrow(dataframe))
dataframe = dataframe[c(3,4,1,2)]

#### Plot NMDS and save the OTU- and tax-tables ####
# Color Palette
col_vector = brewer.pal(n = 12, name = 'Paired')
col_vector[3] = col_vector[4]
col_vector[4] = col_vector[8]
col_vector = col_vector[1:4]
shapes = c(21,22,23,24,25)

# Plot variables
file_name = "NMDS_Hellinger_allOTUs_combined_noSCHF.pdf"
title = "NMDS Plot with the Hellinger transformation for all the OTUs without SchF"
# Summarized plot with all the taxgroups
pdf(file = file_name, width = 9.2, height = 8.3)
  p = ggplot(data = dataframe, aes(x = NMDS1, y = NMDS2, color = taxgroup, fill = taxgroup, shape = Location))
  p = p + geom_point(size = 3, color = "black", alpha = 0.8)
  p = p + guides(fill = guide_legend(override.aes = list(color = col_vector)))
  p = p + scale_fill_manual(values = c("Chlorophyceae" = col_vector[1], 
                                        "Trebouxiophyceae" = col_vector[2], 
                                        "Ulvophyceae" = col_vector[3], 
                                        "Xanthophyceae" = col_vector[4]))
  p = p + scale_shape_manual(values = c("AS14" = shapes[1], 
                                        "AS15" = shapes[2], 
                                        "AM31" = shapes[3], 
                                        "AM06" = shapes[4], 
                                        "AM09" = shapes[5]))
  p = p + ggtitle(title)
  p
dev.off()

# Plot variables
# i {1-4} for the four taxgroups: "Chlorophyceae", "Trebouxiophyceae", "Ulvophyceae" and "Xanthophyceae"
i = 1
names = c("Chlorophyceae", "Trebouxiophyceae", "Ulvophyceae", "Xanthophyceae")
file_name = paste("NMDS_Hellinger_", names[i], "_noSCHF.pdf", sep = "")
title = paste("NMDS Plot with the Hellinger transformation for all the ", names[i], " OTUs without SchF", sep = "")
# Plot for the individual taxgroups
pdf(file = file_name, width = 9.2, height = 8.3)
  p = plist[[i]]
  p = p + aes(fill = taxgroup)
  p = p + scale_color_manual(values = "black")
  p = p + scale_fill_manual(values = c("Chlorophyceae" = col_vector[1], 
                                       "Trebouxiophyceae" = col_vector[2], 
                                       "Ulvophyceae" = col_vector[3], 
                                       "Xanthophyceae" = col_vector[4]))
  p = p + scale_shape_manual(values = c("AS14" = shapes[1],
                                        "AS15" = shapes[2],
                                        "AM31" = shapes[3],
                                        "AM06" = shapes[4],
                                        "AM09" = shapes[5]))
  p = p + ggtitle(title) + theme(plot.title = element_text(size = 12))
  p
dev.off()

# Write the generated OTU- and tax-tables
write.csv(tax_mat, file = "tax_mat.csv")
write.csv(otu_mat, file = "otu_mat.csv")
