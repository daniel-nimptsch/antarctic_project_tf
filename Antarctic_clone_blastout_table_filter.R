################################################################################
# 
# Antarctic_clone_blastout_table_filter.R
#
# Script to delete entries from a blastn_table (BLASTN) depending on a list of
# sequence names provided by the "directions-to-exclude" file.
#
# Written by Daniel Nimptsch
#
################################################################################

#### Packages ####
library(stringr)

#### Working Directory ####
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
setwd(paste(own_cloud_dir, 
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/clones_BLAST/", 
            sep = ""))
list.files()
path = getwd()

# Directions to exclude sequences from the blastout tables
directions = read.csv("directions_to_exclude_sequences_from_fasta.csv", header = TRUE, sep = ";")

# Detect all the blast_tables in the directory
tables = list.files(pattern = "blastn_table")

for (y in 1:length(tables)) {
  # Load the table
  table = read.csv(tables[y], header = TRUE, sep = "\t")
  new_table = matrix(data = NA, ncol = ncol(table))
  # Delete rows from the table corresponding to the directions and create a match_matrix to 
  # verify which items have been deleted
  match_matrix = matrix(data = NA, ncol = 3)
  for (i in 1:nrow(directions)) {
    if (directions[i,2] != "") {
      match = strsplit(strsplit(directions[i,2], "\\:")[[1]][2], "\\|")[[1]][1]
      del_ind = which(table[,1] == match)
      if (is.na(match_matrix[1,1])) {
        match_matrix[1,1] = match
        match_matrix[1,2] = paste(del_ind, collapse = "_")
        match_matrix[1,3] = length(del_ind)
      } else {
        match_matrix = rbind(match_matrix, append(append(match, paste(del_ind, collapse = "_")), length(del_ind)))
      }
      if (length(del_ind) != 0) {
        if (is.na(new_table[1,1])) {
          new_table = table[del_ind,]
        } else {
          new_table = rbind(new_table, table[del_ind,])
        }
        table = table[-del_ind,]
      }
    }
  }
  match_matrix[which(match_matrix[,2] == ""),2] = "no_match"
  colnames(match_matrix) = c("match_string_directions", "row_index_blastout_table", "length_row_index")
  
  # Save the new table
  if (!dir.exists("sorted")) dir.create("sorted")
  name = str_replace(paste(tables[y]), "long_all", "matches")
  name = paste("sorted/", name, sep = "")
  write.table(table, str_replace(name, "matches", "without_matches"), sep = "\t", row.names = FALSE)
  write.table(new_table, str_replace(name, "matches", "only_matches"), sep = "\t", row.names = FALSE)
  write.table(match_matrix, str_replace(name, "blastn_table", "match_record"), sep = "\t", row.names = FALSE)
}

