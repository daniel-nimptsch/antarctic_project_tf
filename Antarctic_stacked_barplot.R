################################################################################
#
# Antarctic_stacked_barplot.R
#
# For the given Table generate a stacked boxplot of the OTU_counts
# for the taxgroups and save it to pdf 
#
# Written by Daniel Nimptsch 
#
################################################################################

#### Packages ####
library(RColorBrewer)
library(ggplot2)
library(stringr)

#### Working Directory ####
path = "data/barplot/"
setwd(path)
list.files()

#### Prepare and tidy data ####
# OTU_counts
for (j in 1:length(list.files(pattern = "all"))) {
  # Load one of the four table for the specific taxgroup
  table = read.csv(file = list.files(pattern = "all")[j], sep = "\t", header = TRUE)
  taxgroup = str_remove(strsplit(list.files(pattern = "all")[j], split = "_")[[1]][3], ".csv")
  table[is.na(table)] = 0
  rownames(table) = table[,1]
  table = table[,-1]
  table = table[,-6]
  colnames(table)[6] = "SchF"
  # Create a temp table only for the taxgroup
  temp_table = matrix(data = NA, nrow = 6, ncol = 3)
  temp_table = as.data.frame(temp_table)
  colnames(temp_table) = c("Location", "Taxgroup", "OTU_count")
  temp_table[,1] = colnames(table)
  temp_table[,2] = taxgroup
  temp_table[,3] = 0
  temp_table[,3] = as.integer(temp_table[,3])
  # Add the OTU_counts
  for (i in 1:nrow(table)) {
    for (y in 1:6) {
      if (table[i,y] > 0) {
        temp_table[y,3] = temp_table[y,3] + 1
      }
    }
  }
  # Add the temp_table to the final_table
  if (exists("final_table")) {
    final_table = rbind(final_table, temp_table)
  } else {
    final_table = temp_table
  } 
}

# Create a final_table_percent for the proportional otu_counts
final_table_percent = final_table
for (i in 1:6) {
  ind = which(final_table[,1] == final_table[i,1])
  sum = sum(final_table[ind,3])
  final_table_percent[ind,3] = final_table_percent[ind,3]/sum
}
final_table_percent[,3] = final_table_percent[,3]*100 

# reads
for (j in 1:length(list.files(pattern = "all"))) {
  # Load one of the four table for the specific taxgroup
  table = read.csv(file = list.files(pattern = "all")[j], sep = "\t", header = TRUE)
  taxgroup = str_remove(strsplit(list.files(pattern = "all")[j], split = "_")[[1]][3], ".csv")
  table[is.na(table)] = 0
  rownames(table) = table[,1]
  table = table[,-1]
  table = table[,-6]
  colnames(table)[6] = "SchF"
  # Create a temp table only for the taxgroup
  reads_table = matrix(data = NA, nrow = 6, ncol = 3)
  reads_table = matrix(data = NA, nrow = 6, ncol = 3)
  reads_table = as.data.frame(reads_table)
  colnames(reads_table) = c("Location", "Taxgroup", "reads")
  reads_table[,1] = colnames(table)
  reads_table[,2] = taxgroup
  reads_table[,3] = 0
  reads_table[,3] = as.integer(reads_table[,3])
  # Add the OTU_counts
  for (i in 1:nrow(table)) {
    for (y in 1:6) {
      if (table[i,y] > 0) {
        reads_table[y,3] = reads_table[y,3] + table[i,y]
      }
    }
  }
  # Add the reads_table to the final_reads
  if (exists("final_reads")) {
    final_reads = rbind(final_reads, reads_table)
  } else {
    final_reads = reads_table
  } 
}

# Create a final_reads_percent for the proportional otu_counts
final_reads_percent = final_reads
for (i in 1:6) {
  ind = which(final_reads[,1] == final_reads[i,1])
  sum = sum(final_reads[ind,3])
  final_reads_percent[ind,3] = final_reads_percent[ind,3]/sum
}
final_reads_percent[,3] = final_reads_percent[,3]*100 


#### Color ####
# Color Palette
col_vector = brewer.pal(n = 12, name = 'Paired')
display.brewer.pal(n = 12, name = 'Paired')
col_vector[3] = col_vector[4]
col_vector[4] = col_vector[8]

#### Plot ####
# Specifications
final_table$Taxgroup <- factor(final_table$Taxgroup, levels = c("Chlorophyceae", "Trebouxiophyceae", "Ulvophyceae", "Xanthophyceae"))
sample_order = c("AM31", "AM09", "AM06", "AS14", "AS15","SchF")

# plot otu counts
pdf(file =  "Boxplot_absolut.pdf", width = 4.2, height = 7.3)
p = ggplot(final_table, aes(x = Location, y = OTU_count))
p = p + geom_bar(aes(color = Taxgroup, fill = Taxgroup), 
                 stat = "identity", position = position_stack())
p = p + geom_text(aes(group = Taxgroup, label = OTU_count),  
                  size = 3, position = position_stack(vjust = 0.5))
p = p + scale_color_manual(values = col_vector) 
p = p + scale_fill_manual(values = col_vector) 
p = p + scale_x_discrete(limits = sample_order)
p = p + ggtitle("Boxplot of the absolut OTU_counts")
p = p + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(size = 11))
p
dev.off()

pdf(file =  "Boxplot_percent.pdf", width = 4.2, height = 7.3)
p = ggplot(final_table_percent, aes(x = Location, y = OTU_count))
p = p + geom_bar(aes(color = Taxgroup, fill = Taxgroup), 
                 stat = "identity", position = position_stack())
p = p + geom_text(aes(group = Taxgroup, label = round(OTU_count, digits = 1)),  
                  size = 3, position = position_stack(vjust = 0.5))
p = p + scale_color_manual(values = col_vector) 
p = p + scale_fill_manual(values = col_vector) 
p = p + scale_x_discrete(limits = sample_order)
p = p + ggtitle("Boxplot of the procentual OTU_counts")
p = p + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(size = 11))
p
dev.off()

# plot reads
final_reads$Taxgroup <- factor(final_reads$Taxgroup, levels = c("Chlorophyceae", "Trebouxiophyceae", "Ulvophyceae", "Xanthophyceae"))
pdf(file =  "Boxplot_reads_absolut.pdf", width = 4.2, height = 7.3)
p = ggplot(final_reads, aes(x = Location, y = reads))
p = p + geom_bar(aes(color = Taxgroup, fill = Taxgroup), 
                 stat = "identity", position = position_stack())
p = p + geom_text(aes(group = Taxgroup, label = reads),  
                  size = 2, position = position_stack(vjust = 0.5))
p = p + scale_color_manual(values = col_vector) 
p = p + scale_fill_manual(values = col_vector) 
p = p + scale_x_discrete(limits = sample_order)
p = p + ggtitle("Boxplot of the absolut reads")
p = p + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(size = 11))
p
dev.off()

pdf(file =  "Boxplot_reads_percent.pdf", width = 4.2, height = 7.3)
p = ggplot(final_reads_percent, aes(x = Location, y = reads))
p = p + geom_bar(aes(color = Taxgroup, fill = Taxgroup), 
                 stat = "identity", position = position_stack())
p = p + geom_text(aes(group = Taxgroup, label = round(reads, digits = 1)),  
                  size = 3, position = position_stack(vjust = 0.5))
p = p + scale_color_manual(values = col_vector) 
p = p + scale_fill_manual(values = col_vector) 
p = p + scale_x_discrete(limits = sample_order)
p = p + ggtitle("Boxplot of the procentual reads")
p = p + theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(size = 11))
p
dev.off()