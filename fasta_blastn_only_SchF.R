################################################################################
#
# fasta_blastn.R
#
# Do a blastn-query with the nt-db for all the fasta-files in the given directory
# and save the blastoutput to a blastoutput-file
#
# Written by Daniel Nimptsch 
#
################################################################################

library(seqinr)

#### Working Directory ####
setwd("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/only_SchF")
path = getwd()
list.files()

# Read the fasta
fasta_file = read.fasta(list.files(pattern = "fasta")[1], as.string = TRUE)

# Read the OTU-Table
otu_mat = read.csv("OtuMatrix_onlySchF.csv", header = TRUE, sep = "\t")

# New fasta
fasta_file_only_SchF = list()
for (i in 1:nrow(otu_mat)) {
  fasta_file_only_SchF = append(fasta_file_only_SchF, fasta_file[which(attributes(fasta_file)$name == otu_mat[i,1])])
}
write.fasta(fasta_file_only_SchF, names = attr(fasta_file_only_SchF, "name"), file.out = "only_SchF.fasta")

# Read the fastas
fastas = list.files(pattern = "fasta")[2]

#### fasta_blastn.sh ####
# Shell script location
script = "/home/pipeline/ownCloud/Arbeit_SAG/Pipeline/r_pipeline_statistics/fasta_blastn.sh"
# iterate through the fastas and do the blastn query
for (i in 1:length(fastas)) {
  path_query = paste(path, "/", fastas[i], sep = "")
  path_out = paste(path, "/", strsplit(fastas[i], split = "\u002Efasta")[[1]][1], ".nt.blastout", sep = "")
  print(paste("BLAST QUERY: ", fastas[i], sep = ""))
  system2("sh", args = c(script, path_query, path_out), wait = TRUE, timeout = 0)
}
