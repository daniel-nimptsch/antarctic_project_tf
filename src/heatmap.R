#
# Antarctic_phyloseq_heatmap.R
#
# For the given tables create phyloseq objects to plot heatmaps of the top 20
# OTUs from the Antarctic project
#
# Written by Daniel Nimptsch
#

# Packages
library(phyloseq)
library(tidyverse)

# Functions
plot_and_save_heatmap <- function(phy) {
    class <- unique(as.data.frame(tax_table(phy)[, 3]))[[1]][1]
    names <- c("C" = "Chlorophyceae", "T" = "Trebouxiophyceae", "U" = "Ulvophyceae", "X" = "Xanthophyceae")
    taxname <- names[which(names(names) == class)]
    sample_order <- c("AM31", "AM09", "AM06", "AS14", "AS15", "SchF")
    taxa_order <- rev(as.vector(taxa_names(phy)))
    title <- paste("Top 20 ", taxname, " taxa", sep = "")

    # Print and save the heatmap to pdf
    # Custom version with edited source code -> y-axis to the right
    p <- plot_heatmap_ypos_right(phy,
        low = "#000033",
        high = "#CCFF66",
        taxa.label = "label_right",
        taxa.order = taxa_order,
        title = title,
        sample.order = sample_order
    )
    p <- p + theme(
        axis.text.x = element_text(size = 10, angle = -45),
        axis.text.y = element_text(size = 10, angle = 0),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 11, hjust = 1)
    )
    p
    ggsave(file.path("./output/heatmaps/", paste(taxname, "_dummy_right.pdf", sep = "")),
        plot = p,
        device = "pdf",
        width = 6.6,
        height = 8.3,
        dpi = 300
    )

    # Default version -> y-axis to the left
    p <- plot_heatmap(phy,
        low = "#000033",
        high = "#CCFF66",
        taxa.label = "label_left",
        taxa.order = taxa_order,
        title = title,
        sample.order = sample_order
    )
    p <- p + theme(
        axis.text.x = element_text(size = 10, angle = -45),
        axis.text.y = element_text(size = 10, angle = 0),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 11, hjust = 1)
    )
    p
    ggsave(file.path("./output/heatmaps/", paste(taxname, "_dummy_left.pdf", sep = "")),
        plot = p,
        device = "pdf",
        width = 6.6,
        height = 8.3,
        dpi = 300
    )
}

# Edit of the plot_heatmap function
source("src/heatmap_custom_functions.R")

# Load data
input_table <- read_delim("data/endversion_2023_Meseta_statistics_1003_OTUs.csv")

taxa_table <- input_table
colnames(taxa_table)[1] <- "otu_id"
taxa_table <- arrange(taxa_table, otu_id)

otu_table <- taxa_table %>%
    select(otu_id, AS14:SchF)
otu_matrix <- otu_table[, -1] %>%
    as.matrix()
rownames(otu_matrix) <- otu_table$otu_id

taxa_matrix <- taxa_table %>%
    filter(!only_SchF) %>%
    select(otu_id, NB, approved_species_taxonomy, class) %>%
    mutate(label_left = str_c(otu_id, NB, approved_species_taxonomy, sep = "_")) %>%
    mutate(label_right = str_c(approved_species_taxonomy, NB, otu_id, sep = "_")) %>%
    select(label_left, label_right, class) %>%
    as.matrix()

rownames(taxa_matrix) <- taxa_table %>%
    filter(!only_SchF) %>%
    select(otu_id) %>%
    as_vector()

phy <- phyloseq(
    otu_table(otu_matrix, taxa_are_rows = T),
    tax_table(taxa_matrix)
)

sample_data(phy) <- sample_data(
    data.frame(
        samples = sample_names(phy),
        sample_type = c("AS14", "AS15", "AM31", "AM06", "AM06", "AM09", "AM09", "SchF"),
        row.names = sample_names(phy)
    )
)

phy_c <- phy %>%
    subset_taxa(class == "C") %>%
    prune_taxa(names(sort(taxa_sums(.), T)[1:20]), .) %>%
    merge_samples("sample_type", mean)

phy_c <- phyloseq(
    otu_table(cbind(otu_table(phy_c), "dummy" = c(rep(0, 5), 4096)), taxa_are_rows = FALSE),
    tax_table(rbind(tax_table(phy_c), "dummy" = c("dummy", "dummy", "dummy"))),
    sample_data(phy_c)
)

phy_t <- phy %>%
    subset_taxa(class == "T") %>%
    prune_taxa(names(sort(taxa_sums(.), T)[1:20]), .) %>%
    merge_samples("sample_type", mean)

phy_u <- phy %>%
    subset_taxa(class == "U") %>%
    prune_taxa(names(sort(taxa_sums(.), T)[1:20]), .) %>%
    merge_samples("sample_type", mean)

phy_x <- phy %>%
    subset_taxa(class == "X") %>%
    prune_taxa(names(sort(taxa_sums(.), T)[1:20]), .) %>%
    merge_samples("sample_type", mean)

plot_and_save_heatmap(phy_c)
plot_and_save_heatmap(phy_t)
plot_and_save_heatmap(phy_u)
plot_and_save_heatmap(phy_x)
