################################################################################
#
# Antarctic_stacked_boxplot.R
#
# For the given Table generate a stacked boxplot for the taxgroups and save it to pdf 
#
# Written by Daniel Nimptsch 
#
################################################################################

#### Packages ####
library(RColorBrewer)
library(ggplot2)

#### Working Directory ####
setwd("/home/pipeline/ownCloud/Arbeit_SAG/Pipeline_Results/Antarctis_1_NGS/Antarctis_1_NGS_2020/only_SchF")
list.files()

table = read.csv(file = list.files(pattern = "Prydzbay")[1], sep = "\t", header = TRUE)
graph_table = table

#### Color ####
# Color Palette
col_vector = brewer.pal(n = 12, name = 'Paired')
col_vector[3] = col_vector[4]
col_vector[4] = col_vector[8]
col_vector[5] = col_vector[12]
black = brewer.pal(n = 9, "Greys")
col_vector[6] = black[9]

#### Plot ####
graph_table$Taxgroup <- factor(graph_table$Taxgroup, levels = 
              c("Chlorophyceae", "Trebouxiophyceae", "Ulvophyceae", "Xanthophyceae", "others", "Streptophyta_Klebsi"))
sample_order = c("PrydzBay", "Meseta_OTUs", "Meseta_SchF_shared_OTUs", "SchF_OTUs")

pdf(file =  "Meseta_SchF_shared_OTUs_PrydzBay.pdf", width = 4.2, height = 8.3)
p = ggplot(graph_table, aes(x = Location, y = OTU_count))
p = p + geom_bar(aes(color = Taxgroup, fill = Taxgroup), 
                 stat = "identity", position = position_stack()) 
p = p + scale_color_manual(values = col_vector) + scale_fill_manual(values = col_vector) + scale_x_discrete(limits = sample_order)
p = p + ggtitle("Meseta_SchF_shared_OTUs_PrydzBay")
p = p + theme(axis.text.x=element_text(angle=30, hjust = 1))
p
dev.off()
