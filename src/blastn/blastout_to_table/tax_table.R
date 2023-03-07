#
# tax_table.R
#
# Written by Daniel Nimptsch
#

get_tax_table_genbank <- function(final_table) {
  
  require(foreach)
  require(doParallel)
  
  message("Get genbank tax_table")
  
  ranks <- c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
  tax_table <- as.data.frame(matrix(data = NA, nrow = nrow(final_table), ncol = 9))
  colnames(tax_table) <- append(ranks, c("seq_ID", "path"), after = 0)
  tax_table$seq_ID <- final_table$seq_ID
  rows <- nrow(tax_table)
  gc()
  registerDoParallel(detectCores()) # use multicore, set to the number of our cores
  tax_table <- foreach(i = 1:rows, .combine = rbind) %dopar% {
    path <- final_table[i, 4]
    taxa <- as.character(unlist(strsplit(path, split = ";")[[1]]))
    if (length(taxa) == 7) {
      tax_table[i, c(3:9)] <- taxa
    } else {
      warning("taxa was not of length 7: ", paste(taxa, collapse = " "))
    }
    tax_table[i, 2] <- path
    return(tax_table[i, ])
  }
  stopImplicitCluster()
  detach("package:doParallel", unload = TRUE)
  detach("package:foreach", unload = TRUE)
  return(tax_table)
}

get_tax_table_silva <- function(final_table) {
  
  get_ranks_path <- function(ind, path, species, tax_table, silva_taxcat) {
    path_length <- length(strsplit(path, ";")[[1]])
    for (y in 1:path_length) {
      path_rank <- strsplit(path, ";")[[1]][1:y]
      path_rank <- paste(path_rank, collapse = ";")
      path_rank <- paste(path_rank, ";", sep = "")
      rank_name <- strsplit(path_rank, ";")[[1]]
      rank_name <- tail(rank_name, n = 1)
      rank_ind <- which(silva_taxcat$V1 == path_rank)
      rank <- silva_taxcat$V3[rank_ind]
      col_rank_ind <- which(colnames(tax_table) == rank)
      tax_table[ind, col_rank_ind] <- rank_name
      tax_table$species[ind] <- species
    }
    return(tax_table)
  }
  
  silva_get_tax_table_foreach <- function(i, final_table, silva_taxmap, silva_taxcat) {
    accession <- strsplit(final_table$accession[i], split = "\\.")[[1]][1]
    silva_taxmap_ind <- which(silva_taxmap$primaryAccession == accession)[1]
    path <- silva_taxmap$path[silva_taxmap_ind]
    species <- silva_taxmap$organism_name[silva_taxmap_ind]
    tax_table <- get_ranks_path(i, path, species, tax_table, silva_taxcat)
    tax_table$path[i] <- paste(path, species, sep = "")
    return(tax_table[i, ])
  }
  
  require(foreach)
  require(doParallel)

  message("Get silva tax_table")

  silva_dir <- Sys.getenv("SILVADB_DIR")
  silva_taxmap <- read.csv(file = file.path(silva_dir, "taxmap_slv_ssu_ref_nr_138.1.txt"), header = TRUE, sep = "\t")
  silva_taxcat <- read.csv(file = file.path(silva_dir, "tax_slv_ssu_138.1.txt"), header = FALSE, sep = "\t")

  ranks <- c(
    "domain",
    "major_clade",
    "superkingdom",
    "kingdom",
    "subkingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "infraphylum",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "superorder",
    "order",
    "suborder",
    "superfamily",
    "family",
    "subfamily",
    "genus",
    "species"
  )

  tax_table <- as.data.frame(matrix(data = NA, nrow = nrow(final_table), ncol = 23))
  colnames(tax_table) <- append(ranks, c("seq_ID", "path"), after = 0)
  tax_table$seq_ID <- final_table$seq_ID
  rows <- nrow(tax_table)

  gc()
  registerDoParallel(detectCores()) # use multicore, set to the number of our cores
  tax_table <- foreach(i = 1:rows, .combine = rbind) %dopar% silva_get_tax_table_foreach(i, final_table, silva_taxmap, silva_taxcat)
  stopImplicitCluster()

  detach("package:doParallel", unload = TRUE)
  detach("package:foreach", unload = TRUE)

  return(tax_table)
}
