################################################################################
# 
# ITSx.R
#
# 
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
            "/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/clones_BLAST/clones_without_matches_ITSx_BLASTN", 
            sep = ""))
list.files()
path = getwd()

######################################################################################################
#### ITSx ####

project_name = list.files(pattern = "fasta")
project_name = strsplit(project_name, split = "\\.")[[1]][1]

itsx = paste(own_cloud_dir, "/programms_daniel/Pipeline/r_pipeline_statistics/itsx.sh", sep = "")
if (!dir.exists("ITSx")) dir.create("ITSx")
itsx_output_header = paste("[ITSx]_", project_name, sep = "")
asv_fasta_file = paste(path, "/", list.files(pattern = "fasta"), sep = "")
system2("bash", args = c(itsx, asv_fasta_file, itsx_output_header), wait = TRUE, timeout = 0)
