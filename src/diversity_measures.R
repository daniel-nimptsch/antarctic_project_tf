#
# Antarctic_phyloseq_diversity_measures.R
#
# For the given tables create phyloseq objects to plot different diversity
# indices
#
# Written by Daniel Nimptsch
#

# Packages
library(RColorBrewer)
library(phyloseq)
library(tidyverse)

# Functions
get_richness <- function(input_table, tax_class) {
  otu_mat <- input_table %>%
    select(OTU_ID, class, AS14:SchF) %>%
    filter(class == tax_class)
  tax_mat <- input_table %>%
    select(OTU_ID, class) %>%
    filter(class == tax_class)
  phy <- phyloseq(
    otu_table(otu_mat[, -c(1, 2)], taxa_are_rows = TRUE),
    tax_table(tax_mat[, -1])
  )
  taxa_names(phy) <- otu_mat$OTU_ID
  richness <- estimate_richness(phy, split = TRUE, measures = c("Shannon", "InvSimpson", "Observed"))
  richness <- bind_cols("sample" = rownames(richness),
                        richness,
                        "taxgroup" = rep(names(taxgroups[taxgroups == tax_class]), nrow(richness))) %>%
    pivot_longer(!c(taxgroup, sample), names_to = "index", values_to = "value")
  return(richness)
}

# Source custom functions
source("pipeline_statistics_custom_phyloseq_functions.R")

# Load data
input_table <- read_delim("./data/2023_Meseta_statistics_1003_OTUs.csv")
taxgroups <- c("Chlorophyceae" = "C", "Trebouxiophyceae" = "T", "Ulvophyceae" = "U", "Xanthophyceae" = "X")

div_table <- tibble()
for (class in taxgroups) {
  div_table = bind_rows(div_table, get_richness(input_table, class))
}

# Color
# Color Palette
col_vector <- brewer.pal(n = 12, name = "Paired")
# display.brewer.pal(n = 12, name = "Paired")
col_vector[3] <- col_vector[4]
col_vector[4] <- col_vector[8]

# Plots
# Sample order for the x-axis
sample_order <- c("AM31", "AM09", "AM914", "AM06", "AM614", "AS14", "AS15", "SchF")

# Save
# Diversity Plots
div_table$index <- factor(div_table$index, levels = c("Observed", "Shannon", "InvSimpson"))
p <- ggplot(div_table, aes(x = sample, y = value)) +
  facet_wrap(~index, scales = "free") + 
  geom_boxplot(alpha = 0.7, coef = 1.5) +
  scale_color_manual(values = col_vector) +
  geom_point(aes(color = taxgroup), size = 3) + 
  theme_bw() + 
  scale_x_discrete(labels = c("AM31-13", "AM09-13", "AM09-14", "AM06-13", "AM06-14", "AS14-14", "AS15-14", "SchF")) + 
  ggtitle("Alpha Diversity Plots for the different taxgroups") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12))
p$data$sample <- factor(p$data$sample, levels = sample_order)

ggsave(
  "./output/alpha_diversity_indices/all_taxgroups_diversity_plots_merged_themeBW.pdf", 
  p, 
  width = 10.2, 
  height = 6.3
)

# Significance
measures <- c("Observed", "Shannon", "InvSimpson")
for (i in 1:length(measures)) {
  measure <- measures[i]
  div_sig <- div_table[which(div_table$index == measure), ]
  div_sig <- as_tibble(div_sig)

  qqnorm(div_sig$value)
  qqline(div_sig$value)
  bartlett.test(value ~ sample, div_sig)

  kruskal_test <- kruskal.test(value ~ sample, div_sig)
  dunn_test = FSA::dunnTest(value ~ sample, div_sig, method = "bh")

  sink(file = str_glue("./output/alpha_diversity_indices/Significance_kruskal_{measure}.txt"))
  print(kruskal_test)
  print("")
  print(dunn_test)
  sink()
}

