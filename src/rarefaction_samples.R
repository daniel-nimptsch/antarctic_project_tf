#
# Antarctic_rarefaction_samples.R
#
#
# Written by Daniel Nimptsch
#

# Packages
library(RColorBrewer)
library(vegan)
library(tidyverse)

# Color
# Color Palette
col_vector <- brewer.pal(n = 8, name = "Paired")

# Rarefaction curve
input_table <- read_delim("data/endversion_2023_Meseta_statistics_1003_OTUs.csv") %>%
    select(AS14:SchF)
legend <- colnames(input_table)
input_table <- input_table %>%
    t() %>%
    as.data.frame() %>%
    mutate_all(function(x) as.numeric(x))
raremax <- min(rowSums(input_table))
filename <- file.path("./output/rarefaction/", "Rarefaction_curve_samples.pdf")

# Plot
pdf(file = filename, width = 11, height = 8)
par(mar = c(5.1, 5.1, 4.1, 9), xpd = FALSE)
p <- rarecurve(input_table, step = 20, sample = raremax, col = col_vector, lwd = 1.5, ylab = "OTUs", label = FALSE)
par(xpd = TRUE)
title("Rarefaction curve of the OTU reads from the samples")
legend("right", inset = c(-0.16, 0), legend = legend, col = col_vector, cex = 0.8, lwd = 2, title = "samples:")
dev.off()

# Alternative
filename <- file.path("./output/rarefaction/", "Rarefaction_curve_samples_alternative.pdf")
pdf(file = filename, width = 11, height = 8)
par(mar = c(5.1, 5.1, 4.1, 9), xpd = FALSE)
p <- rarecurve(input_table, step = 20, col = col_vector, lwd = 1.5, ylab = "OTUs", label = FALSE)
par(xpd = TRUE)
title("Rarefaction curve of the OTU reads from the samples")
legend("right", inset = c(-0.16, 0), legend = legend, col = col_vector, cex = 0.8, lwd = 2, title = "samples:")
dev.off()
