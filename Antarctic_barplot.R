################################################################################
# 
# Antarctic_barplot.R
#
# This Scripts takes the final OTU- and taxonomy-tables and plots barplots of the OTUs
# corresponding to the four taxgroups "Chlorophyceae", "Trebouxiophyceae", "Ulvophyceae"
# and "Xanthophyceae"
#
# Written by Daniel Nimptsch
#
################################################################################

#### Packages ####
library(ggplot2)
library(phyloseq)
library(RColorBrewer)

#### Functions ####
# Create a phyloseq object with a given taxonomy table and a OTU/ASV-table (as matrix)
create_physeq_obj = function(otu_mat, tax_mat) {
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  physeq = phyloseq(OTU, TAX)
  return(physeq)
}

# Transform sample counts to relative abundance in percent
transform_otu_mat_percent = function(otu_mat) {
  for (i in 1:length(otu_mat[1,])) {
    otu_mat[,i] = otu_mat[,i]/colSums(otu_mat)[i]
  }
  otu_mat = otu_mat*100
  return(otu_mat)
}

#### Working Machine ####
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")

#### Working directory ####
setwd(paste(own_cloud_dir, 
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/final_tables_for_grafics/NDMS"
            , sep = ""))
path = getwd()
list.files()

# Read tables
otu_mat = as.data.frame(read.table("otu_mat.csv", sep = ",", row.names = 1, header = TRUE))
otu_mat = type.convert(otu_mat)
tax_mat = as.matrix(read.table("tax_mat.csv", sep = ",", row.names = 1, header = TRUE))

# Exclude the otus not corresponding to the four main taxgroups
otu_mat = otu_mat[-which(tax_mat[,2] == "other"),]
tax_mat = tax_mat[-which(tax_mat[,2] == "other"),]

#### Color ####
# Color Palette
col_vector = brewer.pal(n = 12, name = 'Paired')
col_vector[3] = col_vector[4]
col_vector[4] = col_vector[8]
col_vector = col_vector[1:4]

sample_order = c("AM31", "AM09", "AM06", "AS14", "AS15")
physeq = create_physeq_obj(otu_mat, tax_mat)
pdf(file = "Barplot_taxgroup.pdf", width = 5.2, height = 8.3)
  p = plot_bar(physeq, fill = "taxgroup")
  p$data$taxgroup = factor(p$data$taxgroup,  levels = rev(names))
  p = p + geom_bar(aes_string(color = "taxgroup", fill = "taxgroup"), stat = "identity", position = "stack")
  p = p + scale_fill_manual(values = rev(col_vector), aesthetics = c("fill", "colour"), na.value = "#999999")
  p = p + scale_x_discrete(limits = sample_order)
  p = p + ggtitle("Barplot of the OTU-reads from the Taxgroups")
  p = p + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  p
dev.off()

otu_mat_percent = transform_otu_mat_percent(otu_mat)
physeq = create_physeq_obj(otu_mat_percent, tax_mat)
pdf(file = "Barplot_taxgroup_percent.pdf", width = 5.2, height = 8.3)
  p = plot_bar(physeq, fill = "taxgroup")
  p$data$taxgroup = factor(p$data$taxgroup,  levels = rev(names))
  p = p + geom_bar(aes_string(color = "taxgroup", fill = "taxgroup"), stat = "identity", position = "stack")
  p = p + scale_fill_manual(values = rev(col_vector), aesthetics = c("fill", "colour"), na.value = "#999999")
  p = p + scale_x_discrete(limits = sample_order)
  p = p + ggtitle("Barplot of the percentual OTU-reads from the Taxgroups")
  p = p + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  p
dev.off()