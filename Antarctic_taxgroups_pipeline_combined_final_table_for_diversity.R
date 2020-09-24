################################################################################
# 
# taxgroups_combined_final_table.R
#
# Use the final table from the pipeline and combine it with taxgroups tables
# from Thomas to generate an advanced table with detailed information.
#
# Written by Daniel Nimptsch
#
################################################################################

library(readr)

#### Main ####
#### Working Directory ####
setwd("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/BioDivIndices")
list.files()

#### Selection ####
# Names from the taxgroup tables
taxgroups = list.files(pattern = "final_allOTUs_processed")
# taxgroups: chloro 1, taboux 2, Ulvo 3, Xantho 4

# Load the final table
final_table = read.csv("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/ITS_OLD_LARS/fNMDS_Antarktis_1_data.csv", sep = "\t", header = TRUE)
final_table[is.na(final_table)] = 0

# Load the taxgroup tables
for (t in 1:length(taxgroups)) {
  taxgroup_table = read.csv(taxgroups[t], sep = ",", header = TRUE)
  taxgroup_table = taxgroup_table[,-c(6,7,8,9,10,11)]
  taxgroup_table = cbind(taxgroup_table, matrix(data = 0, nrow = nrow(taxgroup_table), ncol = 8))
  colnames(taxgroup_table)[6:13] = c("AM31", "AM06", "AM0614", "AM09", "AM0914", "AS14", "AS15", "SchF")
  comment(taxgroup_table) = strsplit(strsplit(taxgroups[t], split = "_")[[1]][4], split = "\u002Ecsv")[[1]][1]
  for (i in 1:length(taxgroup_table[,1])) {
    # fill with the information from the final_table
    final_table_ind = which(final_table[,1] == taxgroup_table[i,1])
    grep_pattern = c("31", "6_",  "614", "9_",  "914", "X14", "X15", "SchF")
    sum_vector = c()
    for (y in 1:length(grep_pattern)) {
      sum_vector[y] = sum(final_table[final_table_ind,grep(grep_pattern[y], colnames(final_table))])
    }
    sum_vector[is.na(sum_vector)] = 0
    taxgroup_table[i,6] = sum_vector[1]
    taxgroup_table[i,7] = sum_vector[2]
    taxgroup_table[i,8] = sum_vector[3]
    taxgroup_table[i,9] = sum_vector[4]
    taxgroup_table[i,10] = sum_vector[5]
    taxgroup_table[i,11] = sum_vector[6]
    taxgroup_table[i,12] = sum_vector[7]
    taxgroup_table[i,13] = sum_vector[8]
  }
  file_name = paste("final_allOTUs_for_diversity_", comment(taxgroup_table), ".csv", sep = "")
  write_delim(as.data.frame(taxgroup_table), path = file_name, delim = "\t")
}
