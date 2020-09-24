################################################################################
#
# Antarctic_special_dotplot_bitscore.R
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
setwd("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/taxa_dotplot/dotplot_summary")
list.files()
file = list.files(pattern = "for_dotplot")[2]

# Load the table for the dotplot
table = read.csv(file = file, sep = "\t", header = TRUE)
colnames(table) = c("Taxonomic_name", "Normalized_bitscore", "OTU_ID", "Taxgroup")
table$Normalized_bitscore = as.numeric(gsub(",", "\\.", table$Normalized_bitscore))

#### Color ####
# Color Palette
col_vector = brewer.pal(n = 12, name = 'Paired')
col_vector[3] = col_vector[4]
col_vector[4] = col_vector[8]
col_vector = col_vector[1:4]

# Save the plot to pdf
breaks = seq(0, 1.9, 0.1)
title = "Dotplot of the normalized bitscores from the different genera"

# Dotplot
pdf(file = "Selection_summary_dotplot_jitter.pdf", width = 9.2, height = 7.3)
p = ggplot(table, aes(x = Taxonomic_name, y = Normalized_bitscore))
p = p + facet_grid(cols = vars(Taxgroup), space = "free", scales = "free")
p = p + scale_y_continuous(limits = c(0, 1.85), breaks = seq(0, 1.9, 0.1), minor_breaks = FALSE)
p = p + geom_jitter(width = 0.1, shape = 21, colour = "black", stroke = 0.2, aes(fill = Taxgroup)) + guides(fill = FALSE)
p = p + scale_fill_manual(values = col_vector)
p = p + theme(axis.text.x=element_text(angle = 40, hjust = 1, size = 7.5))
p = p + theme(axis.title.y = element_text(vjust = 2.5))
p = p + ggtitle(title)
p
dev.off()

# Dotplot + Boxplot
pdf(file = "Selection_summary_dotplot_and_boxplot_jitter.pdf", width = 9.2, height = 7.3)
p = ggplot(table, aes(x = Taxonomic_name, y = Normalized_bitscore))
p = p + facet_grid(cols = vars(Taxgroup), space = "free", scales = "free")
p = p + geom_boxplot(alpha = 0.4, lwd = 0.2, outlier.alpha = 0)
p = p + scale_y_continuous(limits = c(0, 1.85), breaks = seq(0, 1.9, 0.1), minor_breaks = FALSE)
p = p + geom_jitter(width = 0.1, shape = 21, colour = "black", stroke = 0.2, aes(fill = Taxgroup)) + guides(fill = FALSE)
p = p + scale_fill_manual(values = col_vector)
p = p + theme(axis.text.x=element_text(angle = 40, hjust = 1, size = 7.5))
p = p + theme(axis.title.y = element_text(vjust = 2.5))
p = p + ggtitle(title)
p
dev.off()
