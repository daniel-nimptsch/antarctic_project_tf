################################################################################
#
# fasta_megablast.R
#
# Do a blastn-query with the nt-db for all the fasta-files in the given directory
# and save the blastoutput to a blastoutput-file
#
# Written by Daniel Nimptsch 
#
################################################################################

#### Working Directory ####
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
setwd(paste(own_cloud_dir, 
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/top_OTUs_allclasses", 
            sep = ""))
list.files()
path = getwd()

# Read the fasta
fastas = list.files(pattern = "top_OTUs_allclasses.fasta")

#### fasta_megablast.sh ####
# Shell script location
script = paste(own_cloud_dir,
               "/programms_daniel/Pipeline/r_pipeline_statistics/fasta_megablast.sh", 
               sep = "")
# iterate through the fastas and do the blastn query
for (i in 1:length(fastas)) {
  fasta_name = strsplit(fastas[i], split = "\u002Efasta")[[1]][1]
  if (length(strsplit(fasta_name, "]")[[1]]) > 1) {
    fasta_name = paste("[megablast]", strsplit(fasta_name, "]")[[1]][2], sep = "")
  } else {
    fasta_name = paste("[megablast]_", fasta_name, sep = "")
  }
  path_query = paste(path, "/", fastas[i], sep = "")
  path_out = paste(path, "/", fasta_name, ".nt.blastout", sep = "")
  print(paste("BLAST QUERY: ", fastas[i], sep = ""))
  system2("sh", args = c(script, path_query, path_out), wait = TRUE, timeout = 0)
}
