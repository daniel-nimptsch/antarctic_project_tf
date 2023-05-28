#
# Antarctic_stacked_barplot.R
#
# For the given Table generate a stacked barplot of the OTU_counts
# for the taxgroups and save it to pdf
#
# Written by Daniel Nimptsch
#

# Packages
library(RColorBrewer)
library(tidyverse)
library(phyloseq)

source("src/phy_to_tidy_table.R")
source("src/table_manipulation.R")

input_table <- read_delim("data/endversion_2023_Meseta_statistics_1003_OTUs.csv")
otu_tab <- input_table %>%
    select(AS14:SchF) %>%
    mutate(meseta = rowMeans(select(., AS14:AM09)))
tax_tab <- input_table %>%
    select(Fig4_boxplot_genus, approved_species_taxonomy, class) %>%
    mutate(path = paste(Fig4_boxplot_genus, approved_species_taxonomy, sep = ";")) %>%
    mutate(class = case_when(
        class == "C" ~ "Chlorophyceae",
        class == "X" ~ "Xanthophyceae",
        class == "T" ~ "Trebouxiophyceae",
        class == "U" ~ "Ulvophyceae",
        .default = "others"
    )) %>%
    select(class, path)

sample_data <- colnames(otu_tab)
sample_data_category <- sample_data
sample_data_category[sample_data_category == "AM614"] <- "AM06"
sample_data_category[sample_data_category == "AM914"] <- "AM09"
sample_data_tab <- bind_cols(sample = sample_data, sample_categories = sample_data_category)
sample_data <- sample_data(sample_data_tab)
sample_names(sample_data) <- sample_data_tab$sample

# Create phyloseq
phy <- phyloseq(
    otu_table(otu_tab, taxa_are_rows = TRUE),
    tax_table(as.matrix(tax_tab)),
    sample_data
)

# Add taxa names to phy
taxa_names(phy) <- input_table$OTU_ID
taxgroup_name <- "class"
tidy_table_list <- get_tidy_table_list(phy, taxgroup_name, "data/")

# Color
# Color Palette
col_vector <- brewer.pal(n = 12, name = "Paired")
col_vector[3] <- col_vector[2]
col_vector[2] <- col_vector[4]
col_vector[4] <- col_vector[8]
col_vector[5] <- "#8b8b8e"
scales::show_col(col_vector)

# Plot
# Specifications
factor_plot_table <- function(table) {
    table$class <- factor(table$class,
        levels = c("Chlorophyceae", "Ulvophyceae", "Trebouxiophyceae", "Xanthophyceae", "others")
    )
    return(table)
}
sample_order <- c("meseta", "AM31", "AM09", "AM914", "AM06", "AM614", "AS14", "AS15", "SchF")

p <- ggplot(
    factor_plot_table(tidy_table_list$tidy_table_absolut_percent),
    aes(x = sample, y = value)
) +
    geom_bar(
        aes(color = class, fill = class),
        stat = "identity", position = position_stack()
    ) +
    geom_text(
        aes(group = class, label = round(value, digits = 0)),
        size = 3, position = position_stack(vjust = 0.5)
    ) +
    scale_color_manual(values = col_vector) +
    scale_fill_manual(values = col_vector) +
    scale_x_discrete(limits = sample_order) +
    ggtitle("Barplot of the relative OTU_counts") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(size = 9))
ggsave("output/stacked_barplot/RelAbsolut.pdf", p, width = 5.2, height = 4.3)

p <- ggplot(
    factor_plot_table(tidy_table_list$tidy_table_percent),
    aes(x = sample, y = value)
) +
    geom_bar(
        aes(color = class, fill = class),
        stat = "identity", position = position_stack()
    ) +
    geom_text(
        aes(group = class, label = round(value, digits = 0)),
        size = 3, position = position_stack(vjust = 0.5)
    ) +
    scale_color_manual(values = col_vector) +
    scale_fill_manual(values = col_vector) +
    scale_x_discrete(limits = sample_order) +
    ggtitle("Barplot of the relative OTU reads") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(size = 9))
ggsave("output/stacked_barplot/RelReads.pdf", p, width = 5.2, height = 4.3)
