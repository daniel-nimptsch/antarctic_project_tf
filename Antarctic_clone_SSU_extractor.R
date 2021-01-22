################################################################################
# 
# clones_SSU_extractor.R
#
# Script used to extract the position from the SSU of the long_clones seq
# and to generate a SSU-fasta from said seq
#
# Written by Daniel Nimptsch
#
################################################################################

#### Packages ####
library(stringr)
library(seqinr)

#### Working Directory ####
own_cloud_dir = Sys.getenv("OWNCLOUD_DIR")
setwd(paste(own_cloud_dir, 
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/clones_BLAST/clones_oncemore", 
            sep = ""))
list.files()
path = getwd()

fasta = list.files(pattern = "fasta")
fasta = read.fasta(fasta)

table = list.files(pattern = "positions")
table = read.csv(table, sep = "\t", header = FALSE)

header = c("seq_name", "length", "SSU", "ITS1", "5.8S", "ITS2", "LSU", "comment")
colnames(table) = header

new_table = matrix(data = NA, nrow = nrow(table), ncol = 3)
colnames(new_table) = c("seq_name", "begin", "end")
new_table[,1] = table$seq_name

for (i in 1:nrow(table)) {
  vector = strsplit(str_remove(table$SSU[i], "SSU: "), "-")[[1]]
  if (length(vector) == 1) {
    new_table[i,c(2,3)] = NA
  } else {
    new_table[i,c(2,3)] = vector
  }
}

for (i in 1:nrow(new_table)) {
  fasta_ind = which(attr(fasta, "name") == new_table[i,1])
  if (is.na(new_table[i,2])) {
    fasta = fasta[-fasta_ind]
  } else {
    fasta_length = length(fasta[[fasta_ind]])
    ssu_begin = as.integer(new_table[i,2])
    ssu_end = as.integer(new_table[i,3])
    fasta[[fasta_ind]] = fasta[[fasta_ind]][ssu_begin:ssu_end] 
  }
}

write.fasta(fasta, attr(fasta, "name"), file.out = "all_clones_SSU.fasta")
