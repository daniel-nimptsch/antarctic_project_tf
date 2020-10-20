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

library(seqinr)

#### Functions ####
# Fill the fasta_mat with information from the fasta file
fasta_fill = function(fasta_mat, fasta_file) {
  fasta_length = getLength(fasta_file)
  for (i in 1:nrow(fasta_mat)) {
    fasta_mat[i,1] = attr(fasta_file[i], "name")
    fasta_mat[i,2] = fasta_length[i]
    fasta_mat[i,3] = fasta_file[[i]][1]
  }
  return(fasta_mat)
}

#### Working Directory ####
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
setwd(paste(own_cloud_dir, 
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/final_OTUs_NGS_Ant1", 
            sep = ""))
list.files()

# Load the final table
final_table = read.csv(paste(own_cloud_dir,"/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/ITS_OLD_LARS/fNMDS_Antarktis_1_data.csv", sep = ""), 
                       sep = "\t", header = TRUE)
final_table[is.na(final_table)] = 0

# Names from the taxgroup tables
taxgroups = list.files(pattern = "final_allOTUs")

# Load th fasta  with the sequences from all the OTUs and generate a table
fasta_file = read.fasta(list.files(pattern = "fasta")[1], as.string = TRUE)
nr_otu = length(fasta_file)
fasta_table = matrix(data = NA, nrow = nr_otu, ncol = 3)
colnames(fasta_table) = c("OTU_ID", "length", "seq")
fasta_table = fasta_fill(fasta_table, fasta_file)

# Load the taxgroup tables
for (t in 1:length(taxgroups)) {
  taxgroup_table = read.csv(taxgroups[t], sep = ";", header = TRUE)
  duplicates = which(duplicated(taxgroup_table[,1]))
  if (length(duplicates) == 0) {
    print("No duplicates were found")
  } else {
    taxgroup_table = taxgroup_table[-duplicates,]
    print("The duplicates were removed:")
    print(duplicates)
  }
  comment(taxgroup_table) = strsplit(strsplit(taxgroups[t], split = "_")[[1]][3], split = "\u002Ecsv")[[1]][1]
  for (i in 1:length(taxgroup_table[,1])) {
    # fill with the information from the fasta
    fasta_ind = which(fasta_table[,1] == taxgroup_table[i,1])
    taxgroup_table[i,4] = fasta_table[fasta_ind,2]
    # fill with the information from the final_table
    final_table_ind = which(final_table[,1] == taxgroup_table[i,1])
    grep_pattern = c("31", "6", "9", "X14", "X15", "SchF")
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
    taxgroup_table[i,5] = sum(sum_vector)
  }
  # write.csv(taxgroup_table, file = paste("processed/final_allOTUs_processed_", comment(taxgroup_table), ".csv", sep = ""), sep = "\t", row.names = FALSE)
}

# Generate a fasta file with the missing OTUs from the Xan taxgroup
taxgroup_table = read.csv(taxgroups[4], sep = ";", header = TRUE)
missing_xantho_fasta = list()
for (i in 1:length(taxgroup_table[,1])) {
  if (taxgroup_table[i,2] == "") {
    otu_fasta_ind = which(names(fasta_file) == taxgroup_table[i,1])
    name_otu_fasta_ind = attr(fasta_file[otu_fasta_ind], "name")
    missing_xantho_fasta[name_otu_fasta_ind] = fasta_file[otu_fasta_ind]
  }
}
# write.fasta(missing_xantho_fasta, names(missing_xantho_fasta), "missing_otu_xantho.fasta")

# Use the information from the BLAST from final_table_missing_otu_xantho.csv to complete the final_allOTUs_Xanthophyceae
blast_table_xantho = read.csv("final_table_missing_otu_xantho.csv", sep = "\t", header = TRUE)
final_table_xantho = read.csv("processed/final_allOTUs_processed_Xanthophyceae.csv", sep = ",", header = TRUE)
for (i in 1:length(final_table_xantho[,1])) {
  if (final_table_xantho[i,2] == "") {
    blast_table_xantho_ind = which(blast_table_xantho[,1] == final_table_xantho[i,1])[1]
    final_table_xantho[i,2] = round(blast_table_xantho[blast_table_xantho_ind,3], digits = 2)
  }
}
# write.csv(final_table_xantho, "processed/final_allOTUs_processed_Xanthophyceae.csv", sep = "\t", row.names = FALSE)

# change the commas in the bitscore for dots
setwd(paste(own_cloud_dir, "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/final_OTUs_NGS_Ant1/processed", 
            sep = ""))
tables = list.files(pattern = ".csv")
for (i in 1:length(tables)) {
  taxgroup_table = read.csv(tables[i], sep = ",", header = TRUE)
  taxgroup_table[,2] = gsub(",", ".", taxgroup_table[,2])
  # write.csv(taxgroup_table, tables[i], sep = "\t", row.names = FALSE)
}

# Load final_table_others_blast
setwd(paste(own_cloud_dir, "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/others_BLAST", sep = ""))
final_table_others = read.csv("final_table_others_BLAST.csv", header = TRUE, sep = "\t")
taxgroup_table = matrix(data = NA, nrow = nrow(fasta_table), ncol = 11)
colnames(taxgroup_table) = c("OTU_ID", "corrected_score", "taxonomy", "seq_length", "total_reads", "AM31", "AM06", "AM09", "AS14", "AS15", "SchF")
taxgroup_table[,1] = fasta_table[,1]
taxgroup_table[,4] = fasta_table[,2]
for (i in 1:length(taxgroup_table[,1])) {
  ind = which(final_table_others[,1] == taxgroup_table[i,1])[1]
  taxgroup_table[i,2] = final_table_others[ind,3]
  taxgroup_table[i,3] = final_table_others[ind,4]
  
  # fill with the information from the final_table
  final_table_ind = which(final_table[,1] == taxgroup_table[i,1])
  grep_pattern = c("31", "6", "9", "X14", "X15", "SchF")
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
  taxgroup_table[i,5] = sum(sum_vector)
}
# write.csv(taxgroup_table, file = "final_others_OTU_table.csv", sep = "\t", row.names = FALSE)
