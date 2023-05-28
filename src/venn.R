#
# Antarctic_Venn.R
#
# Written by Daniel Nimptsch
#

# Packages
library(VennDiagram)
library(tidyverse)
library(eulerr)

# Data
input_table <- read_delim("data/endversion_2023_Meseta_statistics_1003_OTUs.csv")

venn_table <- input_table %>%
    filter(!(class %in% c("Z_unclear", "Z_bryophyte", "Z_fungus", "Z"))) %>%
    mutate(meseta = rowMeans(select(., AS14:AM09)))

otus <- venn_table$OTU_ID
meseta <- otus[which(venn_table$meseta > 0)]
schf <- otus[which(venn_table$SchF > 0)]
list <- list("Meseta" = meseta, "SchF" = schf)

width <- 6
height <- 6

p <- plot(euler(list))
ggsave("output/venn/Venn-Diagramm_euler.pdf", p, width = width, height = height)

# Plot
p <- venn.diagram(list,
    fill = c("lightblue", "green"), alpha = c(0.5, 0.5), lwd = 0,
    filename = NULL,
    print.mode = c("raw", "percent"), fontfamily = "sans", sub.fontfamily = "sans",
    main.fontfamily = "sans", cat.fontfamily = "sans"
)

ggsave("output/venn/Venn-Diagramm.pdf", p, width = width, height = height)
p <- venn.diagram(list,
    fill = c("lightblue", "green"), alpha = c(0.5, 0.5), lwd = 0,
    filename = NULL,
    print.mode = "raw", fontfamily = "sans", sub.fontfamily = "sans",
    main.fontfamily = "sans", cat.fontfamily = "sans"
)
ggsave("output/venn/Venn-Diagramm_no_percent.pdf", p, width = width, height = height)

p <- venn.diagram(list,
    fill = c("darkgreen", "green"), alpha = c(0.5, 0.5), lwd = 0,
    filename = NULL,
    print.mode = c("raw", "percent"), fontfamily = "sans", sub.fontfamily = "sans",
    main.fontfamily = "sans", cat.fontfamily = "sans"
)
ggsave("output/venn/Venn-Diagramm_alternative_color.pdf", p, width = width, height = height)

p <- venn.diagram(list,
    fill = c("darkgreen", "green"), alpha = c(0.5, 0.5), lwd = 0,
    filename = NULL,
    print.mode = "raw", fontfamily = "sans", sub.fontfamily = "sans",
    main.fontfamily = "sans", cat.fontfamily = "sans"
)
ggsave("output/venn/Venn-Diagramm_alternative_color_no_percent.pdf", p, width = width, height = height)
