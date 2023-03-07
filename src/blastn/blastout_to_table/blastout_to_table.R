#
# blastout_to_table.R
#
# For a given blastout(csv) and a corresponding fasta-file generate a
# blast_table with processed information.
#
# Written by Daniel Nimptsch
#

# Function to apply the script to a file
# file and fasta must have the same name structure
# Blast query can be blastn and megablast
blastout_to_table <- function(blastout_file) {

  # Libraries
  require(tidyverse)
  require(seqinr)
  require(R.utils)

  # Functions
  # Source every file in the blastout_to_table path
  source_blastout_to_table_functions <- function() {
    source <- list.files("src/blastn/blastout_to_table/")
    source <- file.path("src/blastn/blastout_to_table", source)
    for (i in 1:length(source)) {
      source(source[i])
    }
  }
  source_blastout_to_table_functions()

  # blast_name and project_name
  # path is relative to project path
  path <- dirname(blastout_file)
  blastout_file <- basename(blastout_file)
  # the blast_name is the second string delimited by "_"
  # assuming blastout_file has naming scheme xxx_blastname_*.nt.blastout
  blast_name <- str_split(blastout_file, "_")[[1]][2]
  # assuming prevous naming scheme the
  # project name should be the * in xxx_blastname_*.nt.blastout
  project_name <- get_project_name(blastout_file)

  # Give the Blastout table colnames
  column_names_blastout <- c(
    "sci_names",
    "bitscore",
    "evalue",
    "coverage",
    "identity",
    "accession",
    "tax_id",
    "subject_title"
  )

  # Load blastout table
  blastouttable <- read.csv(file.path(path, blastout_file), sep = "\t", header = FALSE, col.names = column_names_blastout, quote = "")

  # Load the fasta file
  fasta_file <- list.files(path = path, pattern = paste(project_name, ".fas", sep = ""))
  fasta_file <- read.fasta(file.path(path, fasta_file), as.string = TRUE)

  # Generate a fasta matrix
  # with columns seq_ID, sequence length, and sequence
  message("Get fasta matrix")
  fasta_mat <- get_fasta_mat(fasta_file)
  
  # Get the number of seqs that were blasted
  nr_seq <- length(fasta_file)
  
  # Generate a matrix with the indices of the seq_IDs in the blastout table
  table_indices_mat <- get_table_indices(blastouttable, nr_seq)

  # final_table
  message("Generate final table")
  final_table <- generate_final_table(table_indices_mat, fasta_mat, blastouttable)

  # Filter the final table to only retain the top 10 hits
  final_table <- filter_final_table_top_hits(final_table, nr_seq, fasta_mat, 10)

  # Fixes for different blastdb:
  # Silva
  if (blast_name == "silva") {
    final_table_lineage <- silva_taxa_lineage(final_table)
    # Unite
  } else if (blast_name == "unite") {
    final_table_lineage <- unite_taxa_lineage(final_table)
    # Genbank
  } else if (blast_name == "genbank") {
    final_table_lineage <- genbank_taxa_lineage(final_table)
  }

  # Only_first_hit
  final_table_first_hit <- only_fist_hit(final_table_lineage)

  # Determine the consensus taxonomy
  final_table_consensus <- get_consensus_taxonomy(final_table_lineage, path)

  # Tax Tables
  if (blast_name == "silva") {
    tax_table_first_hit <- get_tax_table_silva(final_table_first_hit)
    tax_table_consensus <- get_tax_table_silva(final_table_consensus)
  } else if (blast_name == "genbank") {
    tax_table_first_hit <- get_tax_table_genbank(final_table_first_hit)
    tax_table_consensus <- get_tax_table_genbank(final_table_consensus)
  }
  
  if (!dir.exists(file.path(path, "rds"))) {
    dir.create(file.path(path, "rds"))
  }
  
  # Save the final_table
  write.table(final_table_first_hit,
    path_table(blast_name, "final-table", "first-hit", "csv", project_name, path),
    sep = ",", row.names = FALSE
  )
  
  saveRDS(final_table_first_hit, path_table(blast_name, "final-table", "first-hit", "rds", project_name, path))

  write.table(final_table_consensus,
    path_table(blast_name, "final-table", "consensus", "csv", project_name, path),
    sep = ",", row.names = FALSE
  )
  
  saveRDS(final_table_first_hit, path_table(blast_name, "final-table", "consensus", "rds", project_name, path))
  
  # Save the final_table_long
  write.table(final_table_lineage,
    path_table(blast_name, "final-table", "long", "csv", project_name, path),
    sep = ",", row.names = FALSE
  )
  
  write.table(tax_table_first_hit,
    path_table(blast_name, "tax-table", "consensus", "csv", project_name, path),
    sep = ",", row.names = FALSE
  )
  
  saveRDS(final_table_first_hit, path_table(blast_name, "tax-table", "consensus", "rds", project_name, path))

  write.table(tax_table_consensus,
    path_table(blast_name, "tax-table", "first-hit", "csv", project_name, path),
    sep = ",", row.names = FALSE
  )
  
  saveRDS(final_table_first_hit, path_table(blast_name, "tax-table", "first-hit", "rds", project_name, path))

  detach("package:seqinr", unload = TRUE)
  detach("package:R.utils", unload = TRUE)
}
