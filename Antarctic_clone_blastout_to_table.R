######################################################################################################
#
# Antarctic_clone_blastout_to_table.R
#
# Script to make a csv-table out of a BLASTN blastoutput file specific to the clones 
# alignment.
#
# Written by Daniel Nimptsch 
#
######################################################################################################

#### Packages ####
library(seqinr)
library(stringr)

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

# Get the start and stop row indices for the up to 100 hits of the table
get_table_indices = function(table, nr_otu) {
  table_indices_mat = matrix(data = NA, nrow = nr_otu, ncol = 3)
  y = 1
  for (i in 1:nrow(table)) {
    if (grepl("*# Query*", table[i,1])) {
      table_indices_mat[y,1] = table[i,1] 
      table_indices_mat[y,2] = i
      y = y + 1
    }
  }
  for (i in 1:nrow(table_indices_mat)) {
    start = as.integer(table_indices_mat[i,2])
    start = start + 4
    table_indices_mat[i,2] = start
    if (i == nr_otu) {
      stop = nrow(table) - 1
      table_indices_mat[i,3] = stop
    } else {
      stop = i + 1 
      stop = as.integer(table_indices_mat[stop,2])
      stop = stop - 2
      table_indices_mat[i,3] = stop
    }
  }
  return(table_indices_mat)
}

# Extract the top hits (hits with top score) and enter them to the final_table and 
# also enter the information from the fasta_mat
generate_final_table = function(table_indices_mat, fasta_mat, table) {
  for (i in 1:nrow(table_indices_mat)) {
    OTU = fasta_mat[i]
    start = as.integer(table_indices_mat[i,2])
    stop = as.integer(table_indices_mat[i,3])
    
    if (start > stop) {
      query_otu_mat = matrix(data = NA, nrow = 1, ncol = 8)
      placeholder_otu_mat = matrix(data = NA, nrow = 1, ncol = 1)
      placeholder_otu_mat[,1] = OTU
      colnames(placeholder_otu_mat) = "OTU_ID"
      colnames(query_otu_mat) = colnames(table)
      query_otu_mat = cbind(placeholder_otu_mat, query_otu_mat)
      query_otu_mat[1,2] = "no_hit"
      placeholder_fasta_mat = matrix(data = NA, nrow = 1, ncol = 3)
      placeholder_fasta_mat[1,1] = fasta_mat[i,1]
      placeholder_fasta_mat[1,2] = fasta_mat[i,2]
      colnames(placeholder_fasta_mat) = c("OTU_ID", "seq_length", "normalized_bitscore")
      query_otu_mat = cbind(placeholder_fasta_mat, query_otu_mat)
    }
    
    else {
      score = table[start, 2]
      y = as.integer(start)
      while ((y < nrow(table)) & y <= stop) {
        y = y + 1
      }
      y = y - 1
      row = y - start + 1
      end = as.integer(y)
      query_otu_mat = matrix(data = NA, nrow = row, ncol = 8)
      placeholder_otu_mat = matrix(data = NA, nrow = row, ncol = 1)
      placeholder_otu_mat[,1] = OTU
      colnames(placeholder_otu_mat) = "OTU_ID"
      query_otu_mat = table[start:end,]
      colnames(query_otu_mat) = colnames(table)
      query_otu_mat = cbind(placeholder_otu_mat, query_otu_mat)
      placeholder_fasta_mat = matrix(data = NA, nrow = nrow(query_otu_mat), ncol = 3)
      for (e in 1:nrow(query_otu_mat)) {
        placeholder_fasta_mat[e,1] = fasta_mat[i,1]
        placeholder_fasta_mat[e,2] = fasta_mat[i,2]
      }
      colnames(placeholder_fasta_mat) = c("OTU_ID", "seq_length", "normalized_bitscore")
      query_otu_mat = cbind(placeholder_fasta_mat, query_otu_mat)
    }
    if (i == 1) {
      final_table = query_otu_mat
    } else {
      final_table = rbind(final_table, query_otu_mat)
    }
  }
  final_table[,3] = as.integer(final_table[,6]) / as.integer(final_table[,2])
  final_table = final_table[,-4]
  return(final_table)
}

#--------------------------------------------------
#### Working Directory ####

# # For all clones:
# own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
# setwd(paste(own_cloud_dir, 
#             "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/clones_BLAST/clones_alignment",
#             sep = ""))
# list.files()
# 
# algal_clones_unmatched:
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
setwd(paste(own_cloud_dir,
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/clones_BLAST/clones_comparison_for_submission/alignment_unmatched",
            sep = ""))
list.files()
#
# matched_otu_clones:
# own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
# setwd(paste(own_cloud_dir,
#             "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/clones_BLAST/clones_comparison_for_submission/alignment_matched",
#             sep = ""))
# list.files()

# Variables
# project_name = "matched_otu_clones"
# file = "matched_otu_clones_megablast_blastoutput"

project_name = "assigned_algal_clones_unmatched"
file = "assigned_algal_clones_unmatched_megablast_blastoutput"
column_names_blastout = c("bitscore", "evalue", "coverage", "identity", "subject_title")

#### Load the Blastout table ####
# Read the blastoutput
table = read.csv(file, sep = "\t", header = FALSE, col.names = column_names_blastout)
# "sci_names",  "bitscore", "evalue", "coverage", "identity", "accession", "tax_id", "subject_title")

# Add columns to the blastoutput so I can reuse the code from the blastoutput_to_table script
table = cbind(rep(NA, nrow(table)), table[,c(1:4)], rep(NA, nrow(table)),rep(NA, nrow(table)), table[,5])
colnames(table) = c("sci_names",  "bitscore", "evalue", "coverage", "identity", "accession", "tax_id", "subject_title")
table[,1] = table[,2]

# Read the corresponding fasta file
fasta_file = read.fasta(list.files(pattern = paste(project_name, ".fasta", sep = ""))[1], as.string = TRUE)
nr_otu = length(fasta_file)
fasta_mat = matrix(data = NA, nrow = nr_otu, ncol = 3)
colnames(fasta_mat) = c("OTU_ID", "length", "seq")
fasta_mat = fasta_fill(fasta_mat, fasta_file)

####  Generate a matrix with the indices ####
table_indices_mat = get_table_indices(table, nr_otu)

#### final_table ####
final_table = generate_final_table(table_indices_mat, fasta_mat, table)

# Now sort the final table:
final_table[,4] = rep(NA, nrow(final_table))
deletion_vector = c()
for (i in 1:nrow(final_table)) {
  # If the seq matched with itself add it to the deletion vector
  # if (final_table[i,1] == final_table[i,11]) {
  #   if (length(deletion_vector) == 0) {
  #     deletion_vector[1] = i
  #   } else {
  #     deletion_vector = append(deletion_vector, i)
  #   }
  # If the matched seq has a nb value below 1.7 add it to the deletion vector  
  # } else if (final_table[i,3] < 1.7) {
  if (final_table[i,3] < 1.75) {
    if (length(deletion_vector) == 0) {
      deletion_vector[1] = i
    } else {
      deletion_vector = append(deletion_vector, i)
    }
  }
}

# Now remove the seq from the deletion vector from the final table
final_table = final_table[-deletion_vector,]

# Just output the unique seqs matched with the library
length(unique(final_table$OTU_ID))

# Remove empty columns
final_table = final_table[,-c(4,9,10)]

# Save the final table to csv
blast_name = "megablast"

file_descriptor = paste("/[", blast_name, "_table]_", sep = "")
write.table(final_table, 
            paste(getwd(), file_descriptor, project_name, ".csv", sep = ""), 
            sep = "\t", row.names = FALSE)
