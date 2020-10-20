################################################################################
# 
# Antarctic_clone_fasta_filter.R
#
# Script to delete entries from a fasta file depending on a list of
# sequence names provided by the "directions-to-exclude" file.
#
# Written by Daniel Nimptsch
#
################################################################################

#### Packages ####
library(stringr)
library(seqinr)

#### Working Directory ####
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
setwd(paste(own_cloud_dir, 
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/clones_BLAST/", 
            sep = ""))
list.files()
path = getwd()

# Directions to exclude sequences from the blastout tables
directions = read.csv("directions_to_exclude_sequences_from_fasta.csv", header = TRUE, sep = ";")

fasta = read.fasta(file = list.files(pattern = "long_all_clones.fasta"))

match_matrix = matrix(data = NA, ncol = 3)
for (i in 1:nrow(directions)) {
  if (directions[i,2] != "") {
    match = strsplit(strsplit(directions[i,2], "\\:")[[1]][2], "\\|")[[1]][1]
    del_ind = which(attr(fasta, "name") == match)
    if (is.na(match_matrix[1,1])) {
      match_matrix[1,1] = match
      match_matrix[1,2] = paste(del_ind, collapse = "_")
      match_matrix[1,3] = length(del_ind)
    } else {
      match_matrix = rbind(match_matrix, append(append(match, paste(del_ind, collapse = "_")), length(del_ind)))
    }
    if (length(del_ind) != 0) {
      fasta = fasta[-del_ind]
    }
  }
}
match_matrix[which(match_matrix[,2] == ""),2] = "no_match"
colnames(match_matrix) = c("match_string_directions", "row_index_blastout_table", "length_row_index")
write.fasta(sequences =  fasta, names = attr(fasta, "name"), file.out = "long_all_clones_without_matches.fasta")
