################################################################################
# 
# blastout_to_table.R
#
# For a given blastout(csv) and a corresponding fasta-file generate a blast_table with processed information
#
# Written by Daniel Nimptsch
#
################################################################################

library(seqinr)
library(taxonomizr)
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
      while ((score == table[y,2] & y < nrow(table)) & y <= stop) {
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

# Use taxonomizr to add the scientific names to the final table
add_taxa_sci_names = function(final_table) {
  # leave as it is:
  accession_sql = Sys.getenv("ACCESSIONDB")
  for (i in 1:nrow(final_table)) {
    taxaID = final_table[i,10]
    if (!is.na(taxaID)) {
      taxa <- getTaxonomy(taxaID, accession_sql, 
                          desiredTaxa = c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
      final_table[i,4] = paste(taxa, collapse = ";")
    }
  }
  return(final_table)
}

# Transform the final_table and show only the first hit
only_fist_hit = function(final_table) {
  final_table_first_hit = matrix(data = NA, nrow = 0, ncol = ncol(final_table))
  y = 1
  for (i in 1:nrow(final_table)) {
    if (i > y) {
      otu_ind = which(final_table[,1] == final_table[i,1])
      final_table_first_hit = rbind(final_table_first_hit, final_table[otu_ind[1],])
      length_otu_ind = length(otu_ind)
      y = otu_ind[length_otu_ind]
    }
  }
  return(final_table_first_hit)
}

# Silva-fix
silva_fix = function(final_table) {
  silva_dir = Sys.getenv("SILVADB_DIR")
  silva_taxmap = read.csv(file = paste(silva_dir, "/", "taxmap_slv_ssu_ref_nr_138.1.txt", sep = ""), header = TRUE, sep = "\t")
  for (i in 1:nrow(final_table)) {
    accession = strsplit(final_table$accession[i], split = "\\.")[[1]][1]
    silva_taxmap_ind = which(silva_taxmap$primaryAccession == accession)
    final_table$tax_id[i] = silva_taxmap$taxid[silva_taxmap_ind]
    final_table$sci_names[i] = silva_taxmap$organism_name[silva_taxmap_ind]
  }
  return(final_table)
}


# Silva-fix alternative taxa
silva_alternative_taxa = function(final_table) {
  silva_dir = Sys.getenv("SILVADB_DIR")
  silva_taxmap = read.csv(file = paste(silva_dir, "/", "taxmap_slv_ssu_ref_nr_138.1.txt", sep = ""), header = TRUE, sep = "\t")
  for (i in 1:nrow(final_table)) {
    accession = strsplit(final_table$accession[i], split = "\\.")[[1]][1]
    silva_taxmap_ind = which(silva_taxmap$primaryAccession == accession)
    final_table$sci_names[i] = paste(silva_taxmap$path[silva_taxmap_ind], silva_taxmap$organism_name[silva_taxmap_ind], sep = "")
  }
  return(final_table)
}

#### Working Directory ####
# important to change:
file = list.files(pattern = ".nt.blastout")[1]
project_name = str_remove(file, ".nt.blastout")
project_name = str_remove(project_name, "\\[blastn\\]_")
project_name = str_remove(project_name, "\\[megablast\\]_")
project_name = str_remove(project_name, "\\[blastn_silva\\]_")
project_name = str_remove(project_name, "\\[megablast_silva\\]_")
blast_name = str_remove(strsplit(file, "\\]")[[1]][1], "\\[")

# Give the Blastout table colnames
column_names_blastout = c("sci_names",  "bitscore", "evalue", "coverage", "identity", "accession", "tax_id", "subject_title")

#### Load the Blastout table ####
table = read.csv(file, sep = "\t", header = FALSE, col.names = column_names_blastout)

#### Load the fasta file and genearte a matrix ####
fasta_file = read.fasta(list.files(pattern = paste(project_name, ".fasta", sep = "")), as.string = TRUE)
nr_otu = length(fasta_file)
fasta_mat = matrix(data = NA, nrow = nr_otu, ncol = 3)
colnames(fasta_mat) = c("OTU_ID", "length", "seq")
fasta_mat = fasta_fill(fasta_mat, fasta_file)

####  Generate a matrix with the indices ####
table_indices_mat = get_table_indices(table, nr_otu)

#### final_table ####
final_table = generate_final_table(table_indices_mat, fasta_mat, table)
# If Blast query: SILVA
if (blast_name == "blastn_silva" || blast_name == "megablast_silva") {
  final_table = silva_fix(final_table)
  final_table_alternative_taxa = silva_alternative_taxa(final_table)
} else {
  final_table_alternative_taxa = add_taxa_sci_names(final_table)
}

#### only_first_hit ####
final_table_first_hit = only_fist_hit(final_table)
final_table_alternative_taxa_first_hit = only_fist_hit(final_table_alternative_taxa)

#### Save the final_table ####
file_descriptor = paste("/[", blast_name, "_table]_", sep = "")
write.table(final_table, 
            paste(getwd(), file_descriptor, project_name, ".csv", sep = ""), 
            sep = "\t", row.names = FALSE)
write.table(final_table_alternative_taxa , 
            paste(getwd(), file_descriptor, "_alternative_taxa_", project_name, ".csv", sep = ""), 
            sep = "\t", row.names = FALSE)
write.table(final_table_first_hit, 
            paste(getwd(), file_descriptor, "_first_hit_", project_name, ".csv", sep = ""), 
            sep = "\t", row.names = FALSE)
write.table(final_table_alternative_taxa_first_hit , 
            paste(getwd(), file_descriptor, "_alternative_taxa_first_hit_", project_name, ".csv", sep = ""), 
            sep = "\t", row.names = FALSE)
