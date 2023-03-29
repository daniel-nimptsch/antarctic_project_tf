#
# taxa_lineage.R
#
# Written by Daniel Nimptsch
#

# Use taxonomizr to add the taxonomic lineage to the final table
# The final_table with all the hits may be required in order to
# perform the function get_consensus_taxonomy()

genbank_taxa_lineage <- function(final_table) {
  # Library
  require(taxonomizr)
  require(foreach)
  require(doParallel)
  
  # Function to get the lineage, recursive
  # it divides the task into chunks of max. 20000 lines
  cycle_get_taxonomy <- function(accession_sql, final_table, start = 1, end = 20000) {
    # Row number required for the loop length
    rows <- nrow(final_table)
    # Making sure the end is at max
    # the rownumber of the final table
    if (end > rows) {
      end <- rows
    }
    # Garbage collector
    gc()
    # Use multicore, set to the number of cores
    registerDoParallel(detectCores())
    final_table[start:end, ] <- foreach(i = start:end, .combine = rbind) %dopar% {
      # Check if there are already determined lineages through the
      # taxid
      if (start > 1) {
        ind <- which(final_table$tax_id[1:(start - 1)] == final_table$tax_id[i])
      } else {
        ind <- c()
      }
      if (!(length(ind) == 0)) {
        final_table$sci_names[i] <- final_table$sci_names[ind]
      } else {
        # Get the taxaID from the final table
        taxaID <- final_table$tax_id[i]
        if (!is.na(taxaID)) {
          taxa <- getTaxonomy(
            taxaID,
            accession_sql,
            desiredTaxa =
              c(
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species"
              )
          )
          # Assign the determined lineage
          final_table$sci_names[i] <- paste(taxa, collapse = ";")
        }
      }
      return(final_table[i, ])
    }
    stopImplicitCluster()
    # If the end is not the final row form the final table repeat
    # task recursivley
    if (end < rows) {
      final_table <- cycle_get_taxonomy(accession_sql, final_table, end + 1, end + 20001)
      # If the end of the final_table is reached return the final table
    } else {
      return(final_table)
    }
  }
  # Message the current progress
  message("Get hits lineage")
  # Get the accession sql database from taxonomyzr from the
  # env variable
  accession_sql <- Sys.getenv("ACCESSIONDB")
  final_table <- cycle_get_taxonomy(accession_sql, final_table)

  detach("package:doParallel", unload = TRUE)
  detach("package:foreach", unload = TRUE)
  detach("package:taxonomizr", unload = TRUE)
  # Finally return the final table with the lineage in the sci_names column
  return(final_table)
}

# Silva
silva_taxa_lineage <- function(final_table) {
  require(foreach)
  require(doParallel)

  message("Get silva hits lineage")

  silva_dir <- Sys.getenv("SILVADB_DIR")
  silva_taxmap <- read.csv(file = str_glue("{silva_dir}/taxmap_slv_ssu_ref_nr_138.1.txt"), header = TRUE, sep = "\t")

  gc()
  registerDoParallel(detectCores()) # use multicore, set to the number of our cores
  rows <- nrow(final_table)

  final_table <- foreach(i = 1:rows, .combine = rbind) %dopar% {
    if (final_table$sci_names[i] != "no_hit") {
      accession <- strsplit(final_table$accession[i], split = "\\.")[[1]][1]
      silva_taxmap_ind <- which(silva_taxmap$primaryAccession == accession)[1]
      path <- silva_taxmap$path[silva_taxmap_ind]
      species <- silva_taxmap$organism_name[silva_taxmap_ind]
      # Fill final table
      final_table$sci_names[i] <- paste(path, species, sep = "")
      final_table$tax_id[i] <- silva_taxmap$taxid[silva_taxmap_ind]
      return(final_table[i, ])
    }
  }

  stopImplicitCluster()

  detach("package:doParallel", unload = TRUE)
  detach("package:foreach", unload = TRUE)

  return(final_table)
}

# Unite
unite_taxa_lineage <- function(final_table) {
  final_table$sci_names <- map_chr(
    final_table$subject_title,
    function(x) str_split(x, "\\|")[[1]][5]
  )
  final_table$sci_names <- map_chr(
    final_table$sci_names,
    function(x) str_remove_all(x, "([:lower:]__)")
  )
  return(final_table)
}

# Microgreen
microgreen_taxa_lineage <- function(final_table, microgreen_tax_table) {
  microgreen_tax_table
  colnames(microgreen_tax_table) <- c("acession", "taxa")
  for (i in 1:nrow(final_table)) {
    if (final_table$accession[i] != "") {
      ind <- which(microgreen_tax_table$acession == final_table$accession[i])
      final_table$sci_names[i] <- microgreen_tax_table$taxa[ind]
    }
  }
  final_table$tax_id <- NA
  final_table$subject_title <- NA
  return(final_table)
}