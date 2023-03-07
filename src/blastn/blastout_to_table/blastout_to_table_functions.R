#
# general_blastout_to_table_functions.R
#
# Written by Daniel Nimptsch
#

# Blastout manipulation functions ---------------------------------------------------------------------------------

# Get the project name out of the file name
get_project_name <- function(file) {
  project_name <- str_remove(file, ".nt.blastout")
  project_name <- str_remove(project_name, "blastn_")
  project_name <- str_remove(project_name, "megablast_")
  project_name <- str_remove(project_name, "silva_")
  project_name <- str_remove(project_name, "unite_")
  project_name <- str_remove(project_name, "genbank_")
  return(project_name)
}

# Fill the fasta_mat with information from the fasta file
fasta_fill <- function(fasta_mat, fasta_file) {
  # Get the seq length of all the seqs of the fasta
  fasta_length <- getLength(fasta_file)
  for (i in 1:nrow(fasta_mat)) {
    # Place the seq name into the seq_ID column
    fasta_mat[i, 1] <- attr(fasta_file[i], "name")
    # Place the seq length
    fasta_mat[i, 2] <- fasta_length[i]
    # Place the actual seq into the table: fasta_mat
    fasta_mat[i, 3] <- fasta_file[[i]][1]
  }
  return(fasta_mat)
}

get_fasta_mat <- function(fasta_file) {
  # Get the number of seqs that were blasted
  nr_seq <- length(fasta_file)
  fasta_mat <- matrix(data = NA, nrow = nr_seq, ncol = 3)
  colnames(fasta_mat) <- c("seq_ID", "length", "seq")
  fasta_mat <- fasta_fill(fasta_mat, fasta_file)
  fasta_mat <- as_tibble(fasta_mat)
  fasta_mat$length <- as.numeric(fasta_mat$length)
  return(fasta_mat)
}

# Get the start and stop row indices for the up to 100 hits of the table
get_table_indices <- function(table, nr_seq) {
  last_line <- grep("# BLAST processed", table$sci_names)
  if (!is_empty(last_line)) {
    table <- table[-last_line, ]
  }
  # Generate a data_frame where the indices will be stored
  table_indices_mat <- matrix(data = NA, nrow = nr_seq, ncol = 3)
  colnames(table_indices_mat) <- c("seq_id", "begin", "end")
  table_indices_mat <- as_tibble(table_indices_mat)
  # Find the begin index for the blastoutput of the individual seq
  table_indices_mat$begin <- which(map(table$sci_names, function(x) {
    grepl("*Query:*", x)
  }) == TRUE)
  # Retrieve the corresponding name
  table_indices_mat$seq_id <- table$sci_names[table_indices_mat$begin]
  table_indices_mat$seq_id <- str_remove(table_indices_mat$seq_id, "# Query: ")
  # Actual begin of the hits
  table_indices_mat$begin <- table_indices_mat$begin + 4
  # Determine the end of the hits
  for (i in 1:nrow(table_indices_mat)) {
    if (i != nrow(table_indices_mat)) {
      table_indices_mat$end[i] <- table_indices_mat$begin[i + 1] - 6
    } else if (i == nrow(table_indices_mat)) {
      table_indices_mat$end[i] <- nrow(table)
    }
  }
  for (i in 1:nrow(table_indices_mat)) {
    if (table_indices_mat$begin[i] > table_indices_mat$end[i]) {
      table_indices_mat$begin[i] <- NA
      table_indices_mat$end[i] <- NA
    }
  }
  return(table_indices_mat)
}

# Extract the top hits (hits with top score) and enter them to the final_table and
# also enter the information from the fasta_mat
generate_final_table <- function(table_indices_mat, fasta_mat, table) {

  # Expand the final tably by the columns "seq_ID", "seq_length" and "normalized_bitscore"
  final_table <- table
  empty_col_ft <- rep(NA, nrow(final_table))
  final_table <- cbind(empty_col_ft, empty_col_ft, empty_col_ft, final_table)
  colnames(final_table)[c(1, 2, 3)] <- c("seq_ID", "seq_length", "normalized_bitscore")

  # Determine the seq_ids that did not get a hit
  no_hit_ind <- which(is.na(table_indices_mat$begin))

  # Determine the begin and end of each query output
  chunks <- which(map(table$sci_names, function(x) {
    grepl("^#", x)
  }) != TRUE)
  chunks <- split(chunks, cumsum(c(1, diff(chunks) != 1)))

  # For each seq query output add the seq length and the normalized bitscore
  # j is used to track the fields that needs to be omited if the seq_id
  # did not generate a hit.
  # chunks only contains information about the seq_ids that have a hit
  j <- 1
  for (i in 1:length(chunks)) {
    # So if the there are actually seq_ids with no hit
    # then the final table wont be filled with that information
    if (!is_empty(no_hit_ind)) {
      # if the i-chunk is equal to the (first) seq_id with no hit
      if (i >= no_hit_ind[j]) {
        # then skip i by adding the value j to y
        y <- i + j
        # if the final no_hit_ind has been passed then do not increment j
        if (j < length(no_hit_ind)) {
          j <- j + 1
        }
        # if the i-chunk is smaller to the j-no_hit_ind then a correction needs
        # to be done in order to properly assign y (since j has been incremented)
      } else if (i < no_hit_ind[j]) {
        y <- i + j - 1
      }
    } else {
      y <- i
    }
    chunks_ind <- chunks[[i]]
    seq_id <- fasta_mat$seq_ID[y]
    seq_length <- fasta_mat$length[y]
    normalized_bitscore <- map_dbl(table$bitscore[chunks_ind], function(x) {
      x / seq_length
    })
    final_table$seq_ID[chunks_ind] <- seq_id
    final_table$seq_length[chunks_ind] <- seq_length
    final_table$normalized_bitscore[chunks_ind] <- normalized_bitscore
  }

  # Delete the comment rows
  final_table <- final_table[-which(map(table$sci_names, function(x) {
    grepl("^#", x)
  }) == TRUE), ]

  # Reset the rownames counter
  rownames(final_table) <- NULL
  return(final_table)
}

path_table <- function(blast_name, table_type, taxonomy_type, file_type, project_name, path) {
  if (file_type == "rds") {
    path <- file.path(path, "rds")
  }
  file_descriptor <- paste(blast_name, "_", table_type, "_", sep = "")
  file.path(
    path,
    paste(file_descriptor,
      taxonomy_type,
      "_",
      project_name,
      ".",
      file_type,
      sep = ""
    )
  )
}
