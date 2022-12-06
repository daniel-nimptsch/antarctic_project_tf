################################################################################
# 
# Antarctic_Venn.R
#
#
# Written by Daniel Nimptsch
#
################################################################################

#### Packages ####
library(VennDiagram)
library(ggVennDiagram)
library(ggplot2)

#### Working directory ####
setwd("data/venn/")
path = getwd()

#### Data ####
otus = paste("otu", 1:963, sep = "")
meseta = otus[1:824]
schf = otus[1:75]
schf = append(schf, otus[825:963])
x <- list("Meseta" = meseta, "SchF" = schf)

#### Plot ####
venn.diagram(x, fill = c("lightblue", "green"), alpha = c(0.5, 0.5), lwd = 0, 
             filename = "Venn-Diagramm.svg",
             print.mode = c("raw", "percent"), fontfamily = "sans", sub.fontfamily = "sans", 
             main.fontfamily = "sans", cat.fontfamily = "sans")

venn.diagram(x, fill = c("lightblue", "green"), alpha = c(0.5, 0.5), lwd = 0, 
             filename = "Venn-Diagramm_no_percent.svg",
             print.mode = "raw", fontfamily = "sans", sub.fontfamily = "sans", 
             main.fontfamily = "sans", cat.fontfamily = "sans")

venn.diagram(x, fill = c("darkgreen", "green"), alpha = c(0.5, 0.5), lwd = 0, 
             filename = "Venn-Diagramm_alternative_color.svg",
             print.mode = c("raw", "percent"), fontfamily = "sans", sub.fontfamily = "sans", 
             main.fontfamily = "sans", cat.fontfamily = "sans")

venn.diagram(x, fill = c("darkgreen", "green"), alpha = c(0.5, 0.5), lwd = 0, 
             filename = "Venn-Diagramm_no_percent_alternative_color.svg",
             print.mode = "raw", fontfamily = "sans", sub.fontfamily = "sans", 
             main.fontfamily = "sans", cat.fontfamily = "sans")

pdf(file = "Venn-Diagramm_ggplot2.pdf", width = 10.2, height = 8.3)
  p = ggVennDiagram(x, n.sides = 6000, color = "darkgreen") 
  p = p + scale_fill_gradient(low = "lightgreen", high = "darkgreen")
  p
dev.off()

pdf(file = "Venn-Diagramm_ggplot2_no_percent.pdf", width = 10.2, height = 8.3)
  p = ggVennDiagram(x, n.sides = 6000, color = "darkgreen", label = "count") 
  p = p + scale_fill_gradient(low = "lightgreen", high = "darkgreen")
  p
dev.off()

setwd("../../")