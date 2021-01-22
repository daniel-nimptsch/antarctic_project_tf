######################################################################################################
#
# Antarctic_clone_submission_selection.R
#
#
# Written by Daniel Nimptsch 
#
######################################################################################################

#### Packages ####
library(seqinr)
library(stringr)

#### Functions ####

# Funcion to get the seq length for a given seq name
fasta_seq_length = function(name) {
  ind = which(attr(master_fasta, "name") == name)
  if (length(ind) == 0) {
    return(NA)
  } else {
  seq_length = seqinr::getLength(master_fasta[ind])
  return(seq_length)
  }
}

# Function to get the OTU-name for a given seq name
get_otu_name = function(name) {
  ind = which(matched_clone_otu_table[,2] == name)
  if (length(ind) == 0) {
    return(NA)
  } else {
    otu_name = matched_clone_otu_table[ind,1]
    return(otu_name)
  }
}

# Function to get the taxonomy from the master_table for a given otu-name
get_taxonomy_matched = function(otu_name) {
  otu_ind = which(master_table$OTU_ID == otu_name)
  taxonomy = master_table$genus_species[otu_ind]
  return(taxonomy)
}

get_taxonomy_unmatched = function(name) {
  ind = which(unmatched_master_blastn_table$clone_name == name)[1]
  taxonomy = unmatched_master_blastn_table$taxonomy[ind]
  return(taxonomy)
}

get_nb_taxonomy_unmatched = function(name) {
  ind = which(unmatched_master_blastn_table$clone_name == name)[1]
  nb = unmatched_master_blastn_table$normalized_bitscore[ind]
  return(nb)
}

get_taxonomy_master_table = function(name) {
  ind = which(master_blastn_table$OTU_ID == name)
  taxonomy = master_blastn_table$sci_names[ind]
  return(taxonomy)
}

get_nb_taxonomy_master_table = function(name) {
  ind = which(master_blastn_table$OTU_ID == name)
  nb = master_blastn_table$normalized_bitscore[ind]
  return(nb)
}

get_its2_length_unmatched = function(name) {
  ind = which(its2_positions_table_unmatched$seq_name == name)[1]
  its2_length = its2_positions_table_unmatched$ITS2[ind]
  its2_length = str_remove(its2_length, "ITS2: ")
  its2_length = strsplit(its2_length, split = "-")[[1]]
  its2_length = as.integer(its2_length[2]) - as.integer(its2_length[1])
  return(its2_length)
}

# Create a table with the clone seq names and the correspinding OTUs
get_matched_clone_otu_table = function() {
  matched_clone_otu_table = matrix(data = NA, nrow = 0, ncol = 2)
  for (i in 1:nrow(master_table)) {
    otu_name = master_table[i,1]
    clone_names = master_table[i,17]
    if (!is.na(clone_names)) {
      clone_names = strsplit(clone_names, ":")[[1]][-1]
      clone_names = sub("\\|.*", "", clone_names)
      matched_clone_otu_table = rbind(matched_clone_otu_table, cbind(rep(otu_name, length(clone_names)), clone_names))
    }
  }
  colnames(matched_clone_otu_table)[1] = "OTU_ID"
  return(matched_clone_otu_table)
}

# Function to add more columns to the blastn_table
extend_and_fill_blastn_table = function(blastn_table) {
  blastn_table = cbind(blastn_table[,c(1:2)], rep(NA, nrow(blastn_table)), rep(NA, nrow(blastn_table)), blastn_table[,c(3:8)],
                       rep(NA, nrow(blastn_table)), rep(NA, nrow(blastn_table)), rep(NA, nrow(blastn_table)))
  colnames(blastn_table)[c(3,4,11,12,13)] = c("matched_otu", "group", "subject_title_seq_length", "subject_title_matched_otu", "best_candidate")
  for (i in 1:nrow(blastn_table)) {
    blastn_table$matched_otu[i] = get_otu_name(blastn_table$OTU_ID[i])
    blastn_table$subject_title_seq_length[i] = fasta_seq_length(blastn_table$subject_title[i])
    blastn_table$subject_title_matched_otu[i] = get_otu_name(blastn_table$subject_title[i])
  }
  return(blastn_table)
}

# Function to group similar seq in the blastn_table
group_blastn_table = function(blastn_table) {
  for (i in 1:nrow(blastn_table)) {
    name = blastn_table$OTU_ID[i]
    if (i == 1) {
      blastn_table$group[i] = name
    } else {
      group_present = which(blastn_table$group == name)[1]
      if (!is.na(group_present)) {
        blastn_table$group[i] = blastn_table$group[group_present]
      } else if (is.na(group_present)) {
        other_group_present = which(blastn_table$subject_title == name)
        other_group_present = other_group_present[which(!is.na(blastn_table$group[other_group_present]))]
        other_group_present = other_group_present[1]
        if (!is.na(other_group_present)) {
          blastn_table$group[i] = blastn_table$group[other_group_present]
        } else {
          blastn_table$group[i] = name
        }
      }
    }
  }
  return(blastn_table)
}

# Function to determine the best candidate for each individual seq
compare_best_candidate = function(blastn_table) {
  for (i in 1:nrow(blastn_table)) {
    if (blastn_table$seq_length[i] > blastn_table$subject_title_seq_length[i]) {
      blastn_table$best_candidate[i] = blastn_table$OTU_ID[i]
    } else {
      blastn_table$best_candidate[i] = blastn_table$subject_title[i]
    }
  }
  return(blastn_table)
}

# Function to determine the best group candidate
get_best_group_candidate_blast_table = function(blastn_table) {
  temp_blastn_table = blastn_table
  temp_blastn_table = temp_blastn_table[order(temp_blastn_table$group),]
  i = 1
  while (i < nrow(temp_blastn_table)) {
    ind_group = which(temp_blastn_table$group == temp_blastn_table$group[i])
    if (length(ind_group) > 1) {
      max = max(sapply(temp_blastn_table$best_candidate[ind_group], fasta_seq_length))
      best_group_candidate = attr(which(sapply(temp_blastn_table$best_candidate[ind_group], fasta_seq_length) == max)[1], "name")
      temp_blastn_table$best_candidate[ind_group] = best_group_candidate
    }
    i = ind_group[length(ind_group)] + 1  
  }
  temp_blastn_table = temp_blastn_table[order(as.integer(row.names(temp_blastn_table))),]
  blastn_table = temp_blastn_table
  return(blastn_table)
}

# Function to add the columns {"otu_not_present_in_best", "final_best"} to the blastn table and fill them with seq-names
# corresponding to the clone OTUs missing from the best group candidate
add_column_otu_not_present = function(blastn_table_matched) {
  blastn_table_matched = cbind(blastn_table_matched, rep(NA, nrow(blastn_table_matched)), rep(NA, nrow(blastn_table_matched)))
  colnames(blastn_table_matched)[c(length(blastn_table_matched) - 1, length(blastn_table_matched))] = c("otu_not_present_in_best", "final_best")
  blastn_table_matched$final_best = blastn_table_matched$best_candidate
  for (i in 1:length(otu_not_present)) {
    ind_otu = which(blastn_table_matched$matched_otu == otu_not_present[i])
    max = max(sapply(blastn_table_matched$OTU_ID[ind_otu], fasta_seq_length))
    best_group_candidate = attr(which(sapply(blastn_table_matched$OTU_ID[ind_otu], fasta_seq_length) == max)[1], "name")
    blastn_table_matched$otu_not_present_in_best[which(blastn_table_matched$OTU_ID == best_group_candidate)[1]] = best_group_candidate
    blastn_table_matched$final_best[which(blastn_table_matched$OTU_ID == best_group_candidate)[1]] = best_group_candidate
  }
  return(blastn_table_matched)
}

#### Working Directory ####
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
setwd(paste(own_cloud_dir, 
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/clones_BLAST/clones_comparison_for_submission",
            sep = ""))
list.files()

#### Load the relevant tables ####
# Table with the clone seq matched with OTUs
blastn_table_matched = read.csv(list.files(pattern = "\\[megablast_table\\]_matched"), sep = "\t", header = TRUE)
# Table with the clone seq not corresponding to any OTUs
blastn_table_unmatched = read.csv(list.files(pattern = "\\[megablast_table\\]_assigned"), sep = "\t", header = TRUE)
# The final_table with all the information corresponding to the Antarctic OTUs
master_table = read.csv(list.files(pattern = "finMeseta"), sep = ";", header = TRUE)
master_table$clone.name[which(master_table$clone.name == "")] = NA
# All the clone seqs
master_fasta = read.fasta(list.files(pattern = "long_all_clones_formatted.fasta"))
# Clone matched fasta
matched_fasta = read.fasta(list.files(pattern = "matched_otu_clones.fasta"))
# Clone unmatched fasta
unmatched_fasta = read.fasta(list.files(pattern = "assigned_algal_clones_unmatched.fasta"))
# Blastn_table from all the long_clones blast
master_blastn_table = read.csv(list.files(pattern = "\\[blastn_table\\]_alternative"), sep = "\t", header = TRUE)
# unmatched_master_blastn_table
unmatched_master_blastn_table = read.csv("assigned_algal_clones_unmatched.csv", sep = ";", header = FALSE)
unmatched_master_blastn_table = unmatched_master_blastn_table[,c(1:5)] 
colnames(unmatched_master_blastn_table) = c("taxonomy", "clone_name", "seq_length", "normalized_bitscore", "blastn_taxonomy")
unmatched_master_blastn_table$clone_name = sub("\\|.*", "", unmatched_master_blastn_table$clone_name)
# its2_positions_table_unmatched
its2_positions_table_unmatched = read.csv(list.files(pattern = "positions"), sep = "\t", header = TRUE)

#### Submission selection by comparison of seq properties ####
# Matched
# Create a table with the clone seq names and the correspinding OTUs
matched_clone_otu_table = get_matched_clone_otu_table()

# Add more columns to the blastn table
# {"matched_otu", "group", "subject_title_seq_length", "subject_title_matched_otu", "best_candidate"}
# Comparison by: matched_otu & seq length
blastn_table_matched = extend_and_fill_blastn_table(blastn_table_matched)
# Group the similar seqs and add the group name to the group column
blastn_table_matched = group_blastn_table(blastn_table_matched)
# Compare each individual seqs with thier match and decide the best candidate
blastn_table_matched = compare_best_candidate(blastn_table_matched)
# Compare the best candidates from each group and decide a best group candidate
blastn_table_matched = get_best_group_candidate_blast_table(blastn_table_matched)

# Which OTUs a represented with the best candidates
best_candidate_otu = sapply(blastn_table_matched$best_candidate, get_otu_name)
best_candidate_otu = unique(best_candidate_otu)
# Which are all the clone OTUs
all_matched_otu = unique(matched_clone_otu_table[,1])
# Which clone OTUs are missing from the best candidates
otu_not_present =  setdiff(all_matched_otu, best_candidate_otu)

# Add a new column to complete the best candidates with other seqs that a corresponding to the missing OTUs
blastn_table_matched = add_column_otu_not_present(blastn_table_matched)

# Test again if there are some clone OTUs missing
best_candidate_otu = sapply(blastn_table_matched$final_best, get_otu_name)
best_candidate_otu = unique(best_candidate_otu)
otu_not_present =  setdiff(all_matched_otu, best_candidate_otu)

# Unmatched
# Repeat the procedure for the unmatches seq
blastn_table_unmatched = extend_and_fill_blastn_table(blastn_table_unmatched)
blastn_table_unmatched = group_blastn_table(blastn_table_unmatched)
blastn_table_unmatched = compare_best_candidate(blastn_table_unmatched)
blastn_table_unmatched = get_best_group_candidate_blast_table(blastn_table_unmatched)
blastn_table_unmatched = blastn_table_unmatched[,-c(3,12)]

# Extra
# matched
fasta_ind_matched = sapply(unique(blastn_table_matched$final_best), function(x){which(attr(master_fasta, "name") == x)})
unique_seq_matched = master_fasta[fasta_ind_matched]
length(unique(blastn_table_matched$OTU_ID))
length(unique(blastn_table_matched$final_best))

# unmatched
fasta_ind_unmatched = sapply(unique(blastn_table_unmatched$best_candidate), function(x){which(attr(master_fasta, "name") == x)})
unique_seq_unmatched = master_fasta[fasta_ind_unmatched]
length(unique(blastn_table_unmatched$OTU_ID))
length(unique(blastn_table_unmatched$best_candidate))

# unique_seq_table
# matched
unique_seq_table_matched = blastn_table_matched
unique_seq_table_matched = unique_seq_table_matched[which(!duplicated(unique_seq_table_matched$final_best)),]
unique_seq_table_matched = as.data.frame(cbind(unique_seq_table_matched$final_best, unique_seq_table_matched$seq_length, unique_seq_table_matched$matched_otu, rep(NA, nrow(unique_seq_table_matched))))
colnames(unique_seq_table_matched) = c("unique_seq", "seq_length", "matched_otu", "taxonomy")
for (i in 1:nrow(unique_seq_table_matched)) {
  unique_seq_table_matched$seq_length[i] = fasta_seq_length(unique_seq_table_matched$unique_seq[i])
  unique_seq_table_matched$matched_otu[i] = get_otu_name(unique_seq_table_matched$unique_seq[i])
  unique_seq_table_matched$taxonomy[i] =  get_taxonomy_matched(unique_seq_table_matched$matched_otu[i])
}

# unmatched
unique_seq_table_unmatched = blastn_table_unmatched
unique_seq_table_unmatched = unique_seq_table_unmatched[which(!duplicated(unique_seq_table_unmatched$best_candidate)),]
unique_seq_table_unmatched = as.data.frame(cbind(unique_seq_table_unmatched$best_candidate, 
                                                 unique_seq_table_unmatched$seq_length, 
                                                 rep(NA, nrow(unique_seq_table_unmatched)),
                                                 rep(NA, nrow(unique_seq_table_unmatched)),
                                                 rep(NA, nrow(unique_seq_table_unmatched)),
                                                 rep(NA, nrow(unique_seq_table_unmatched)), 
                                                 rep(NA, nrow(unique_seq_table_unmatched))))
colnames(unique_seq_table_unmatched) = c("unique_seq", "seq_length", "ITS2_length", "taxonomy_assigned", "nb_taxonomy_assigned", "blastn_taxonomy", "nb_blastn_taxonomy")
for (i in 1:nrow(unique_seq_table_unmatched)) {
  unique_seq_table_unmatched$seq_length[i] = fasta_seq_length(unique_seq_table_unmatched$unique_seq[i])
  unique_seq_table_unmatched$taxonomy_assigned[i] =  get_taxonomy_unmatched(unique_seq_table_unmatched$unique_seq[i])
  unique_seq_table_unmatched$nb_taxonomy_assigned[i] = get_nb_taxonomy_unmatched(unique_seq_table_unmatched$unique_seq[i])
  unique_seq_table_unmatched$blastn_taxonomy[i] = get_taxonomy_master_table(unique_seq_table_unmatched$unique_seq[i])
  unique_seq_table_unmatched$nb_blastn_taxonomy[i] = get_nb_taxonomy_master_table(unique_seq_table_unmatched$unique_seq[i])
  unique_seq_table_unmatched$ITS2_length[i] = get_its2_length_unmatched(unique_seq_table_unmatched$unique_seq[i])
}

# Test
setdiff(attr(unique_seq_matched, "name"), unique_seq_table_matched$unique_seq)
setdiff(attr(unique_seq_unmatched, "name"), unique_seq_table_unmatched$unique_seq)

#### Save the results ####
# Save both the comparison tables
write.table(blastn_table_matched, "[comparison_table]_matched_clones.csv", sep = "\t", row.names = FALSE)
write.table(blastn_table_unmatched, "[comparison_table]_unmatched_clones.csv", sep = "\t", row.names = FALSE)

# Save both the unique_seq_tables
write.table(unique_seq_table_matched, "[unique_seq_table]_matched.csv", sep = "\t", row.names = FALSE)
write.table(unique_seq_table_unmatched, "[unique_seq_table]_unmatched.csv", sep = "\t", row.names = FALSE)

# Save fastas from the unique seq (uniqe best candidates)
write.fasta(unique_seq_matched, attr(unique_seq_matched, "name"), file.out = "unique_seq_matched.fasta")
write.fasta(unique_seq_unmatched, attr(unique_seq_unmatched, "name"), file.out = "unique_seq_unmatched.fasta")

