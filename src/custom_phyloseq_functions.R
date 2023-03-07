#### Functions ####

# With given path: load and format taxmat (as matrix)
load_tax_mat = function(path_to_tax_mat) {
  tax_mat = read.csv(path_to_tax_mat, header = TRUE, sep = "\t")
  row.names(tax_mat) = tax_mat[,1]
  tax_mat = as.matrix(tax_mat[,2:ncol(tax_mat)])
  return(tax_mat)
}

# With given path: load and format otumat (as matrix)
load_otu_mat = function(path_to_otu_mat) {
  otu_mat = read.csv(path_to_otu_mat, header = TRUE, sep = "\t")
  row.names(otu_mat) = otu_mat[,1]
  otu_mat = as.matrix(otu_mat[,2:ncol(otu_mat)])
  otu_mat[is.na(otu_mat)] = 0
  return(otu_mat)
}

# Transform sample counts to relative abundance in percent
transform_otu_mat_percent = function(otu_mat) {
  for (i in 1:length(otu_mat[1,])) {
    otu_mat[,i] = otu_mat[,i]/colSums(otu_mat)[i]
  }
  otu_mat = otu_mat*100
  return(otu_mat)
}

# Create a phyloseq object with a given taxonomy table and a OTU/ASV-table (as matrix)
create_physeq_obj = function(otu_mat, tax_mat) {
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  physeq = phyloseq(OTU, TAX)
  return(physeq)
}

# Rarefaction custom otu_matrix, spreaded
create_otu_mat_spread = function(otu_mat) {
  otu_mat_spread = as.data.frame(otu_mat)
  otu_mat_spread = base::t(otu_mat_spread)
  otu_mat_spread[is.na(otu_mat_spread)] = 0
  return(otu_mat_spread)
}

# Pruned and transformed phyloseq object for diversity measurements 
create_physeq_obj_prunded_diversity = function(physeq) {
  physeq.prune = physeq
  otu_table(physeq.prune)[is.na(otu_table(physeq.prune))] = 0 
  physeq.prune = prune_species(speciesSums(physeq.prune) > 0, physeq.prune)
  return(physeq.prune)
}

# Phyloseq Top N 
top_n_physeq_obj = function(top_number, physeq) {
  topN <- names(sort(taxa_sums(physeq), decreasing = TRUE))[1:top_number]
  physeq.topN <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
  physeq.topN <- prune_taxa(topN, physeq.topN)
  return(physeq.topN)
}

# Phyloseq Top N Location Normalized
top_n_normalized_physeq_obj = function(physeq.topN) {
  physeq.topN.norm = transform_sample_counts(physeq.topN, function(x) x/sum(x))
  physeq.topN.norm_m = merge_samples(physeq.topN.norm, "Location", fun = mean)
  sample_data(physeq.topN.norm_m)$Location <- factor(sample_names(physeq.topN.norm_m))
  physeq.topN.norm_m = transform_sample_counts(physeq.topN.norm_m, function(x) 100 * x/sum(x))
  return(physeq.topN.norm_m)
}

# custom barplot nocolor
custom_barplot_nocolor = function(phy, rank) {
  p = plot_bar(phy, fill=paste(rank, sep = ""))
  p + geom_bar(aes_string(color = paste(rank), fill=paste(rank)), stat="identity", position="stack")
}

# custom barplot
custom_barplot = function(phy, rank, col_vector) {
  p = plot_bar(phy, fill=paste(rank, sep = ""))
  p + geom_bar(aes_string(color=paste(rank), fill=paste(rank)), stat="identity", position="stack") +
    scale_fill_manual(values=col_vector, aesthetics = c("fill", "colour"), na.value="#999999")
}

# custom barplot sample data nocolor
custom_barplot_sample_data_nocolor = function(phy, rank, sample_data) {
  p = plot_bar(phy, fill=paste(rank, sep = ""), x = paste(sample_data, sep = ""))
  p + geom_bar(aes_string(color=paste(rank), fill=paste(rank)), stat="identity", position="stack") + 
    ylab("Relative abundance")
}

# custom barplot sample data
custom_barplot_sample_data = function(phy, rank, col_vector, sample_data) {
  p = plot_bar(phy, fill=paste(rank, sep = ""), x = paste(sample_data, sep = ""))
  p + geom_bar(aes_string(color=paste(rank, sep = ""), fill=paste(rank, sep = "")), stat="identity", position="stack") + 
    ylab("Relative abundance") +
    scale_fill_manual(values=col_vector, aesthetics = c("fill", "colour"), na.value="#999999")
}