#
# Antarctic_Venn.R
#
# Written by Daniel Nimptsch
#

# Packages
library(VennDiagram)
library(ggVennDiagram)
library(tidyverse)

# Data
input_table <- read_delim("data/endversion_2023_Meseta_statistics_1003_OTUs.csv")

venn_table <- input_table %>%
    filter(!(class %in% c("Z_unclear", "Z_bryophyte", "Z_fungus", "Z"))) %>%
    mutate(meseta = rowMeans(select(., AS14:AM09)))

otus <- venn_table$OTU_ID
meseta <- otus[which(venn_table$meseta > 0)]
schf <- otus[which(venn_table$SchF > 0)]
list <- list("Meseta" = meseta, "SchF" = schf)

# Plot
venn.diagram(list,
    fill = c("lightblue", "green"), alpha = c(0.5, 0.5), lwd = 0,
    filename = "output/venn/Venn-Diagramm.svg",
    print.mode = c("raw", "percent"), fontfamily = "sans", sub.fontfamily = "sans",
    main.fontfamily = "sans", cat.fontfamily = "sans"
)

venn.diagram(list,
    fill = c("lightblue", "green"), alpha = c(0.5, 0.5), lwd = 0,
    filename = "output/venn/Venn-Diagramm_no_percent.svg",
    print.mode = "raw", fontfamily = "sans", sub.fontfamily = "sans",
    main.fontfamily = "sans", cat.fontfamily = "sans"
)

venn.diagram(list,
    fill = c("darkgreen", "green"), alpha = c(0.5, 0.5), lwd = 0,
    filename = "output/venn/Venn-Diagramm_alternative_color.svg",
    print.mode = c("raw", "percent"), fontfamily = "sans", sub.fontfamily = "sans",
    main.fontfamily = "sans", cat.fontfamily = "sans"
)

venn.diagram(list,
    fill = c("darkgreen", "green"), alpha = c(0.5, 0.5), lwd = 0,
    filename = "output/venn/Venn-Diagramm_no_percent_alternative_color.svg",
    print.mode = "raw", fontfamily = "sans", sub.fontfamily = "sans",
    main.fontfamily = "sans", cat.fontfamily = "sans"
)

p <- ggVennDiagram(list, n.sides = 6000, color = "darkgreen")
p <- p + scale_fill_gradient(low = "lightgreen", high = "darkgreen")
ggsave("output/venn/Venn-Diagramm_ggplot2.pdf", p, width = 10.2, height = 8.3)

p <- ggVennDiagram(list, n.sides = 6000, color = "darkgreen", label = "count")
p <- p + scale_fill_gradient(low = "lightgreen", high = "darkgreen")
ggsave("output/venn/Venn-Diagramm_ggplot2_no_percent.pdf", p, width = 10.2, height = 8.3)
