#
# fasta_blastn.R
#
# Do a blastn-query with the nt-db for all the fasta-files in the given directory
# and save the blastoutput to a blastoutput-file
#
# Written by Daniel Nimptsch
#

fasta_blastn <- function(query, db = "genbank") {
  # Blastn
  if (db == "genbank") {
    script <- "src/blastn/blastn_shell_scripts/fasta_blastn.sh"
    decorator <- "blastn_genbank"

    # Microgreen
  } else if (db == "microgreen") {
    script <- "src/blastn/blastn_shell_scripts/fasta_blastn_microgreen.sh"
    decorator <- "blastn_microgreen"

    # Silva
  } else if (db == "silva") {
    script <- "src/blastn/blastn_shell_scripts/fasta_blastn_silva.sh"
    decorator <- "blastn_silva"

    # Unite
  } else if (db == "unite") {
    script <- "src/blastn/blastn_shell_scripts/fasta_blastn_unite.sh"
    decorator <- "blastn_unite"
  }

  # BLASTN
  # Iterate through the fastas and do the blastn query
  for (i in 1:length(query)) {
    file_query <- query[i]
    path <- dirname(file_query)
    fasta_name <- basename(query[i])

    # Set the name of the output to [blastn]_xxx
    fasta_name <- str_remove(fasta_name, "\\.fasta")
    fasta_name <- paste(decorator, "_", fasta_name, sep = "")

    # Create the output directory [blastn]_xxx and copy the the fasta
    dir_out <- file.path(path, fasta_name)
    if (!dir.exists(dir_out)) {
      dir.create(dir_out)
    }
    file.copy(file_query, dir_out)

    # Set the name of the actual blstoutput file
    path_out <- file.path(dir_out, paste(fasta_name, ".nt.blastout", sep = ""))

    # Execute the BLASTN shell script
    system2(
      "sh",
      args = c(script, file_query, path_out),
      wait = TRUE,
      timeout = 0
    )
  }
}
