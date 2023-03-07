source("src/pipeline_scripts_dn.R")
blastn_standalone("data/raw/all_1003_OTUs_Antarctis1_NGS_sequences.fasta", "genbank")
blastout_to_table_standalone("data/processed/all_1003_OTUs_Antarctis1_NGS_sequences/", "genbank")
