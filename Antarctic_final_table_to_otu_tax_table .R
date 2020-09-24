#### Packages ####
library(taxonomizr)
library(readr)

#### Path ####
# path = "/home/soye/Volume/OwnCloud/Arbeit_SAG"
path = "/home/pipeline/ownCloud/Arbeit_SAG"

#### Source custom functions ####
setwd(paste(path, "/Pipeline/r_pipeline_statistics", sep = ""))
source("final_table_to_otu_tax_table_functions.R")

#### Working directory ####
setwd(paste(path, "/Pipeline_Results/Antarctis_1_NGS/ITS_OLD_LARS", sep = ""))
accession_sql = "/home/pipeline/bigdata/accessiondb/accessionTaxa.sql"

#### Pipeline specifications ####
marker = "ITS"

#### Custom functions only for the Antarctis project ####

# Generate an alternative OTU-Table for the taxagroups of the Antarctis-project
antarctis_alternative_otutable = function(otu_mat, final_table, indices) {
  otu_mat[is.na(otu_mat)] = 0
  otu_mat = cbind(otu_mat, matrix(data = NA, ncol = 1, nrow = nrow(otu_mat)))
  for (i in 1:nrow(otu_mat)) {
    otu_mat[i,17] = mean(otu_mat[i,indices])
  }
  otu_mat = otu_mat[,17]
  otu_mat = as.matrix(otu_mat)
  colnames(otu_mat) = "sample_sums"
  tax_cols = matrix(data = NA, ncol = 6, nrow = nrow(otu_mat))
  colnames(tax_cols) = c("Chlorophyceae", "Trebouxiophyceae", "Ulvophyceae", "Xanthophyceae", "unidentified_Chlorophyta", "other")
  rownames(tax_cols) = rownames(otu_mat)
  for (i in 1:nrow(otu_mat)) {
    str_tax = as.vector(strsplit(final_table$taxgroup[i], split = "_")[[1]])
    str_lenght = length(str_tax)
    # if there is only one coding character:
    if(str_tax[2] == "") {
      if(str_tax[1] == "G") {
        tax_cols[i, 5] = otu_mat[i]
      } else if(str_tax[1] == "F") {
        tax_cols[i, 6] = otu_mat[i]
      } else if(str_tax[1] == "I") {
        tax_cols[i, 4] = otu_mat[i]
      } else if(str_tax[1] == "C") {
        tax_cols[i, 1] = otu_mat[i]
      }
      # if there are two coding characters:  
    } else if(str_tax[2] == "U") {
        tax_cols[i, 3] = otu_mat[i]
    } else if(str_tax[2] == "C") {
        tax_cols[i, 1] = otu_mat[i]
    } else if(str_tax[2] == "T") {
        tax_cols[i, 2] = otu_mat[i]
    } else if (str_tax[2] == "X") {
        tax_cols[i, 4] = otu_mat[i]
    } else if(str_tax[1] == "B") {
      tax_cols[i, 6] = otu_mat[i]
    # if the questions did not ask for everything:
    } else {
      print(paste("i:", i, ": ", final_table$taxgroup[i], " strtax1: -", str_tax[1], "-", " strtax2: -", str_tax[2], "-", sep = ""))
      print(length(str_tax))
    }
  }
  return(tax_cols)
}

# Function to count the otu numbers corresponding to the taxroups
create_taxgroup_otu_counts = function(alternative_otu_mat) {
  taxgroup_otu_counts = matrix(data = 0, ncol = 1, nrow = 6)
  colnames(taxgroup_otu_counts) = "Value"
  rownames(taxgroup_otu_counts) = colnames(alternative_otu_mat[,1:6])
  for (i in 1:nrow(alternative_otu_mat)) {
    taxgroup = (which(!is.na(alternative_otu_mat[i,])))
    if(taxgroup == 1) {
      taxgroup_otu_counts[1,1] = taxgroup_otu_counts[1,1] + 1 
    } else if (taxgroup == 2) {
      taxgroup_otu_counts[2,1] = taxgroup_otu_counts[2,1] + 1 
    } else if (taxgroup == 3) {
      taxgroup_otu_counts[3,1] = taxgroup_otu_counts[3,1] + 1 
    } else if (taxgroup == 4) {
      taxgroup_otu_counts[4,1] = taxgroup_otu_counts[4,1] + 1 
    } else if (taxgroup == 5) {
      taxgroup_otu_counts[5,1] = taxgroup_otu_counts[5,1] + 1 
    } else if (taxgroup == 6) {
      taxgroup_otu_counts[6,1] = taxgroup_otu_counts[6,1] + 1 
    }
  }
  Total_no_other = sum(taxgroup_otu_counts[1:5,1])
  taxgroup_otu_counts = rbind(taxgroup_otu_counts, Total_no_other)
  return(taxgroup_otu_counts)
}

# Get the indices for the OTUs of the otu-mat that are only present in SchF
get_indices_SchF = function(otu_mat) {
  y = 1
  otu_mat_indices = c()
  for (i in 1:nrow(otu_mat)) {
    if(all(is.na(otu_mat[i,1:14]))) {
      otu_mat_indices[y] = i
      y = y + 1
    }
  }
  return(otu_mat_indices)
}


#########################################################################################

# MAIN ---------------------------

#### Create the matrices {otu_mat & tax_mat} out of "final_table.csv" ####

# Load file
final_table <- read_tsv("fNMDS_Antarktis_1_data.csv")

# Determine the last column of the samples
end_sample_column = "clones"
end_sample_column = which(colnames(final_table) == end_sample_column) - 1

# Create a OTU-matrix
otu_mat <- as.matrix(final_table[,1:end_sample_column])
row.names(otu_mat) <- otu_mat[,1]
otu_mat <- otu_mat[,2:ncol(otu_mat)]
otu_mat <- apply(otu_mat, c(1,2), as.numeric)
NamesList <- list(rownames(otu_mat), colnames(otu_mat))
otu_mat <- matrix(data = otu_mat, ncol = ncol(otu_mat), nrow = nrow(otu_mat), dimnames = NamesList)

# Remove OTU that are only SchF
otu_mat_indices = get_indices_SchF(otu_mat)
# otu_mat_onlySchF = otu_mat[otu_mat_indices,]
# Remove SchF only OTUs
# otu_mat = otu_mat[-otu_mat_indices,]
# otu_mat[is.na(otu_mat)] = 0

# Keep SchF only OTUs
otu_mat = otu_mat[otu_mat_indices,]
otu_mat[is.na(otu_mat)] = 0
otu_mat = cbind(otu_mat, final_table[otu_mat_indices,29])

which(colnames(final_table) == "taxgroup")
colnames(final_table)[29]
nrow(otu_mat)

# Cols for otu_mat
col_SchF = grep("SchF.", colnames(otu_mat))
col_XAN = grep("*Xan$" , colnames(otu_mat))
col_XAN = col_XAN[-which((col_XAN %in% col_SchF) == TRUE)]
col_ITS4 = grep("*ITS4$" , colnames(otu_mat))
col_ITS4 = col_ITS4[-which((col_ITS4 %in% col_SchF) == TRUE)]

# Alternative OTU-matrix for taxagroups of the Antarctis-project
alternative_otu_mat_ITS = antarctis_alternative_otutable(otu_mat, final_table, col_ITS4)
alternative_otu_mat_XAN = antarctis_alternative_otutable(otu_mat, final_table, col_XAN)
alternative_otu_mat = alternative_otu_mat_ITS
alternative_otu_mat[,4] = alternative_otu_mat_XAN[,4]

# Count OTU for taxgroups
# taxgroup_otu_counts = create_taxgroup_otu_counts(alternative_otu_mat)

# Create an accession-matrix
if(marker == "EUK" || marker == "ITS") {
  acc_mat <- as.matrix(final_table[,c(1,which(colnames(final_table) == "accession"))])
  row.names(acc_mat) <- acc_mat[,1]
  acc_mat <- as.matrix(acc_mat[,2])
  colnames(acc_mat) <- "accession"
  acc_mat[,1] <- sub('\\..*','', acc_mat[,1])
}

# Create the matrix "bac_blast_str"
if (marker == "BAC") {
  bac_blast_str <- as.matrix(final_table[,c(1,which(colnames(final_table) == "taxonomy"))])
  row.names(bac_blast_str) <- bac_blast_str[,1]
  bac_blast_str <- as.matrix(bac_blast_str[,2])
  colnames(bac_blast_str) <- "taxonomy"
}

# Create a taxonomy-matrix
tax_mat <- matrix(data = NA, nrow = nrow(otu_mat), ncol = 8)
rownames(tax_mat) <- rownames(otu_mat)
colnames(tax_mat) <- c("Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
if(marker == "EUK" || marker == "ITS") {
  tax_mat = taxMatGetTaxonomyItsEuk(tax_mat, acc_mat, accession_sql)
} else if (marker == "BAC") {
  tax_mat = taxMatGetTaxonomyBac(tax_mat, bac_blast_str, accession_sql)
}

# Fill tha NAs with taxonomic names from other ranks
tax_mat = taxMatFill(tax_mat)

#### Save everything ####
if(!dir.exists("R_Statistik")) dir.create("R_Statistik")
write.table(otu_mat, "R_Statistik/OtuMatrix_onlySchF.csv", sep = "\t", row.names = TRUE, col.names = NA)
write.table(alternative_otu_mat, "R_Statistik/Alt_mean_new_OtuMatrix.csv", sep = "\t", row.names = TRUE, col.names = NA)
write.table(taxgroup_otu_counts, "R_Statistik/taxgroup_otu_counts.csv", sep = "\t", row.names = TRUE, col.names = NA)
write.table(tax_mat, "R_Statistik/TaxonomyMatrix.csv", sep = "\t", row.names = TRUE, col.names = NA)