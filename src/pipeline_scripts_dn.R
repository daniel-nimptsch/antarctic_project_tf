#
# pipeline_scripts_dn.R
#
#
# Written by Daniel Nimptsch
#

# Libraries -------------------------------------------------------------------------------------------------------

# Library used in all the computations
library(tidyverse)

# Functions -------------------------------------------------------------------------------------------------------

# Full blast standalone
# This function may be used in order to do a blastn search with the specified database
# Possible: genbank, unite and silva.
# The fasta_path is the path of the fasta location. The results will be put
# into the data/processed/ directory.
# The blastout_to_table function is also executed in order to generate
# readable tables and in order to provide the normalized bitscore aswell as
# the consensus taxonomy.
blastn_standalone <- function(fasta_path, db = "genbank") {
  pipeline_scripts_path = getwd()

  # The dataset name is the filename of the fasta
  dataset_name <- basename(fasta_path)
  dataset_name <- str_remove(dataset_name, ".fasta")
  dataset_path <- str_glue("data/processed/{dataset_name}/")

  # Create a folder in data/processed/ with tha dataset_name
  # and copy the fasta into it
  if (!dir.exists(dataset_path)) {
    dir.create(dataset_path)
  }
  file.copy(fasta_path, dataset_path)

  # Fasta blastn
  source("src/blastn/fasta_blastn.R")
  fasta <- basename(fasta_path)
  query <- file.path(getwd(), fasta_path)
  if (!file.exists(query)) stop("Query fasta does not exist!")
  fasta_blastn(query, db)
}

blastout_to_table_standalone <- function(dataset_path, db = "genbank") {
  # For each of the named databases execute
  # Blastout to table
  for (i in 1:length(db)) {
    source("src/blastn/blastout_to_table/blastout_to_table.R")
    blastout_file <- get_blastout_path(dataset_path, db[i])
    if (!file.exists(blastout_file)) stop("Blastoutput does not exist!")
    blastout_to_table(blastout_file)
  }
}