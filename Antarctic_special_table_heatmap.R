################################################################################
#
# Antarctic_special_table_heatmap.R
#
# For the given tables create phyloseq objects to plot heatmaps of the top 20
# OTUs from the Antarctic project
#
# Written by Daniel Nimptsch 
#
################################################################################

library(readr)

#### Working directory ####
setwd("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/heatmaps")
file = list.files(pattern = "special")[2]

# Load the final table
final_table = read.csv("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/ITS_OLD_LARS/fNMDS_Antarktis_1_data.csv", sep = "\t", header = TRUE)
final_table[is.na(final_table)] = 0

special = read.csv(list.files(pattern = "special_heatmap")[1], sep = ";", header = TRUE)
colnames(special)[10] = "SchF"

for (i in 1:nrow(special)) {
  final_table_ind = which(final_table[,1] == special[i,1])
  special[i,10] = sum(final_table[final_table_ind, grep("SchF", colnames(final_table))])
}

write_delim(x = special, path = file, delim = "\t")
