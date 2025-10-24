library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)

# Define methods and corresponding colors
methods <- c("C-SIDE", "Celina", "STANCE","spVC")
colors <- c("#0073C2FF","#EFC000FF","#CD534CFF","#868686FF")
method_colors <- setNames(colors, methods)

my_theme <- theme(
  # Overall plot background
  plot.background   = element_blank(),
  panel.background  = element_blank(),
  panel.grid        = element_blank(),
  panel.border      = element_rect(color = "black", linewidth = 0.5, fill = NA),
  axis.line         = element_blank(),
  
  # Axes
  axis.text         = element_text(color = "black", size = 7),
  axis.ticks        = element_line(color = "black", linewidth = 0.5),
  axis.title        = element_text(color = "black", size = 7),
  axis.title.x      = element_text(size = 7),
  axis.title.y      = element_text(size = 7),
  
  # Legend
  legend.background = element_blank(),
  legend.key        = element_blank(),
  legend.key.size   = unit(0.2, "in"), 
  legend.key.width  = unit(0.2, "in"),     
  legend.box.margin = margin(0, 0, 0, 0),
  legend.margin     = margin(0, 0, 0, 0),
  legend.text       = element_text(color = "black", size = 7),
  legend.text.align = 0,
  legend.title      = element_text(color = "black", size = 7),
  
  # Facet titles
  strip.text        = element_text(color = "black", size = 7),
  
  # Overall plot title
  plot.title        = element_text(color = "black", size = 7)
)
