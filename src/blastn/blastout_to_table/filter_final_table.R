#
# filter_final_table.R
#
# Written by Daniel Nimptsch
#

# Filtering functions ---------------------------------------------------------------------------------------------

# Filter the final table by keeping only the hits with the maximal
# normalized bitscore
filter_final_table_max_normalized_bitscore <- function(final_table, nr_seq, fasta_mat) {

  # Empty vector for indices of the maximal normalized bitscore hits
  vec_ind_max_nor_bit <- c()

  # Go through the seqs
  for (i in 1:nr_seq) {
    ind <- which(final_table$seq_ID == fasta_mat$seq_ID[i])

    # Determine the maximal normalized bitscore of the hits of the given seq
    max_nor_bit <- max(final_table$normalized_bitscore[ind])
    ind_max_nor_bit <- which(final_table$normalized_bitscore[ind] >= max_nor_bit)
    ind_max_nor_bit <- ind[ind_max_nor_bit]
    vec_ind_max_nor_bit <- append(vec_ind_max_nor_bit, ind_max_nor_bit)
  }

  # Return the final table with only the hits with the maximal normalized bitscore value
  final_table <- final_table[vec_ind_max_nor_bit, ]
  return(final_table)
}

# Filter the final table with only the hits that have a identity above 90 percent
filter_final_table_identity <- function(final_table) {
  final_table <- final_table[which(final_table$identity >= 90), ]
  return(final_table)
}

# Filter the final table with a specified number of top hits
filter_final_table_top_hits <- function(final_table, nr_seq, fasta_mat, nr_hits) {
  vec_top_hits <- c()
  for (i in 1:nr_seq) {
    ind <- which(final_table$seq_ID == fasta_mat$seq_ID[i])
    if (!is_empty(ind)) {
      if (length(ind) < nr_hits) {
        vec_top_hits <- append(vec_top_hits, ind)
      } else {
        ind <- ind[1]
        ind <- seq(ind, ind + nr_hits - 1, 1)
        vec_top_hits <- append(vec_top_hits, ind)
      }
    }
  }
  final_table <- final_table[vec_top_hits, ]
  return(final_table)
}

# Transform the final_table and show only the first hit
only_fist_hit <- function(final_table) {
  final_table_first_hit <- distinct(final_table, seq_ID, .keep_all = TRUE)
  return(final_table_first_hit)
}