#
# Antarctic_dotplot_bitscore_summarized_clone.R
#
# Script used to plot the normalized-bitscore values derived from the BLAST query.
#
#
# Written by Daniel Nimptsch
#

# Packages
library(RColorBrewer)
library(tidyverse)

# Functions
plot_dotplot_boxplot <- function(dotplot_table, facet = TRUE, bw = FALSE) {
    p <- ggplot(dotplot_table, aes(x = Taxonomic_name, y = Normalized_bitscore))
    if (facet) p <- p + facet_grid(cols = vars(Taxgroup), space = "free_x", scales = "free_x")
    if (bw) p <- p + theme_bw()
    p <- p + scale_y_continuous(
        limits = limits, breaks = breaks, minor_breaks = FALSE,
        sec.axis = sec_axis(~., breaks = breaks)
    ) +
        geom_boxplot(alpha = 0.4, lwd = 0.2, outlier.alpha = 0) +
        geom_jitter(
            data = dotplot_table[dotplot_table$Clone == FALSE, ], width = 0.1, shape = 21, colour = "black",
            stroke = 0.2, aes(fill = Taxgroup)
        ) +
        guides(fill = FALSE) +
        geom_jitter(
            data = dotplot_table[dotplot_table$Clone == TRUE, ], width = 0.1, shape = 23, colour = "black",
            stroke = 0.2, fill = "yellow"
        ) +
        guides(fill = FALSE) +
        scale_fill_manual(values = col_vector) +
        theme(
            strip.text.x = element_text(size = 6),
            axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
            panel.grid.major = element_line(size = 0.2),
            axis.title.y = element_text(vjust = 2.5),
        ) +
        ggtitle(title) +
        ylab("normalized bitscore (NC)")
    return(p)
}

# Color Palette
get_colors <- function() {
    col_vector <- brewer.pal(n = 12, name = "Paired")
    col_vector[3] <- col_vector[4]
    col_vector[4] <- col_vector[8]
    col_vector <- col_vector[1:4]
    return(col_vector)
}

# Load input table
input_table <- read_delim("./data/2023_Meseta_statistics_1003_OTUs.csv")

dotplot_table <- input_table %>%
    select(Fig4_boxplot_genus, NC, OTU_ID, OTU_with_clone, class) %>%
    mutate(NC = as.numeric(gsub(",", "\\.", NC))) %>%
    mutate(OTU_with_clone = !is.na(OTU_with_clone)) %>%
    mutate(class = replace(class, class == "C", "Chlorophyceae")) %>%
    mutate(class = replace(class, class == "T", "Trebouxiophyceae")) %>%
    mutate(class = replace(class, class == "U", "Ulvophyceae")) %>%
    mutate(class = replace(class, class == "X", "Xanthophyceae")) %>%
    filter(!(input_table$class %in% c("Z_unclear", "Z_bryophyte", "Z_fungus"))) %>%
    arrange(Fig4_boxplot_genus)
colnames(dotplot_table) <- c("Taxonomic_name", "Normalized_bitscore", "OTU_ID", "Clone", "Taxgroup")

# Save the plot to pdf
# breaks for the y-axis
breaks <- seq(0, 1.85, 0.1)
# limits for the y-axis
limits <- c(0.1, 1.85)
# Reorder the x-axis labels to the order specified by the .csv
dotplot_table$Taxonomic_name <- factor(dotplot_table$Taxonomic_name, levels = unique(dotplot_table$Taxonomic_name))

# Color
col_vector <- get_colors()
title <- "Dotplot of the normalized bitscores from the different genera"
p <- plot_dotplot_boxplot(dotplot_table, bw = TRUE)
ggsave("./output/dotplot/dotplot_boxplot.pdf", p, width = 12.2, height = 7.3)
# p <- plot_dotplot_boxplot(dotplot_table)
# ggsave("./output/dotplot/Summary_dotplot_and_boxplot_jitter_with_clones.pdf", p, width = 12.2, height = 7.3)

for (taxgroup in unique(dotplot_table$Taxgroup)) {
    title <- str_glue("Dotplot of the normalized bitscores from {taxgroup}")
    if (taxgroup == "Xanthophyceae") {
        width <- 2.2
        col_vector <- get_colors()[4]
    } else if (taxgroup == "Ulvophyceae") {
        width <- 3.2
        col_vector <- get_colors()[3]
    } else if (taxgroup == "Trebouxiophyceae") {
        width <- 6.2
        col_vector <- get_colors()[2]
    } else if (taxgroup == "Chlorophyceae") {
        width <- 6.2
        col_vector <- get_colors()[1]
    }
    p <- dotplot_table %>%
        filter(Taxgroup == taxgroup) %>%
        plot_dotplot_boxplot(facet = FALSE, bw = TRUE)
    ggsave(str_glue("./output/dotplot/dotplot_boxplot_{taxgroup}.pdf"), p, width = width, height = 4.8)
}

title <- str_glue("Dotplot of the normalized bitscores from Ulvo and Xantho")
col_vector <- get_colors()[3:4]
p <- dotplot_table %>%
    filter(Taxgroup %in% c("Ulvophyceae", "Xanthophyceae")) %>%
    plot_dotplot_boxplot(facet = TRUE, bw = TRUE)
ggsave(str_glue("./output/dotplot/dotplot_boxplot_Ulvo_Xantho.pdf"), p, width = 3.2, height = 4.8)
