################################################################################
#
# Antarctic_dotplot_bitscore.R
#
# 
#
# Written by Daniel Nimptsch 
#
################################################################################

#### Packages ####
library(RColorBrewer)
library(ggplot2)

#### Working Directory ####
setwd("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/taxa_dotplot/dotplots")
list.files()
file = list.files(pattern = "for_dotplot")[4]
name = strsplit(file, split = "_")[[1]][1]

# Load the table for the dotplot
table = read.csv(file = file, sep = ";", header = TRUE)
colnames(table) = c("Taxonomic_name", "Normalized_bitscore", "OTU_ID")
table$Normalized_bitscore = as.numeric(gsub(",", "\\.", table$Normalized_bitscore))

#### Color ####
# Color Palette
col_vector = brewer.pal(n = 12, name = 'Paired')
col_vector[3] = col_vector[4]
col_vector[4] = col_vector[8]

if(name == "Chlorophyceae") {
  col_vector = col_vector[1]
  ratio = 1
  width = 7.2
} else if (name == "Trebouxiophyceae") {
  col_vector = col_vector[2]  
  ratio = 0.7
  width = 9.2
} else if (name == "Ulvophyceae") {
  col_vector = col_vector[3]
  ratio = 0.7
  width = 7.2
  taxonnames = unique(table$Taxonomic_name)
} else if (name == "Xanthophyceae") {
  col_vector = col_vector[4]
  ratio = 0.27
  width = 7.2
  taxonnames = unique(table$Taxonomic_name)
}

# Save the plot to pdf
breaks = seq(0, 1.9, 0.1)
title = paste("Dotplot of the normalized bitscores from the different genera from ", name, sep = "")

# Dotplot
pdf(file = paste(name, "_dotplot_jitter.pdf", sep = ""), width = width, height = 7.3)
p = ggplot(table, aes(x = Taxonomic_name, y = Normalized_bitscore))
p = p + scale_y_continuous(limits = c(0, 1.85), breaks = seq(0, 1.9, 0.1), minor_breaks = FALSE)
p = p + geom_jitter(width = 0.1, shape = 21, colour = "black", fill = col_vector, stroke = 0.2)
p = p + theme(axis.text.x=element_text(angle = 40, hjust = 1))
p = p + theme(axis.title.y = element_text(vjust = 2.5))
p = p + ggtitle(title)
if(name == "Ulvophyceae") {
  p = p + scale_x_discrete(limits = taxonnames)
}
p
dev.off()

# Dotplot + Boxplot
pdf(file = paste(name, "_dotplot_and_boxplot_jitter.pdf", sep = ""), width = width, height = 7.3)
p = ggplot(table, aes(x = Taxonomic_name, y = Normalized_bitscore))
p = p + scale_y_continuous(limits = c(0, 1.85), breaks = seq(0, 1.9, 0.1), minor_breaks = FALSE)
p = p + geom_boxplot(alpha = 0.4, lwd = 0.2, outlier.alpha = 0)
p = p + geom_jitter(width = 0.1, shape = 21, colour = "black", fill = col_vector, stroke = 0.2)
p = p + theme(axis.text.x = element_text(angle = 40, hjust = 1))
p = p + theme(axis.title.y = element_text(vjust = 2.5))
p = p + ggtitle(title)
if(name == "Ulvophyceae") {
  p = p + scale_x_discrete(limits = taxonnames)
}
p
dev.off()
