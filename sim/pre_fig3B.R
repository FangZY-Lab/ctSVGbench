data_info=read.csv('./sim_data_info1.csv')
data_info <- data_info %>% 
  select(
    !resolution     
  )
merged_data <- merge(auc, data_info, by = "dt", all.x = TRUE)

library(ggplot2)
library(ggnewscale)  
library(RColorBrewer) 
library(dplyr)        

merged_data <- merged_data %>%
  mutate(method_pattern = paste(pattern, methods, sep = ": ")) %>% 
  complete(dt, method_pattern) %>%
  group_by(dt) %>%
  fill(resolution, tech_platform, species, tissue, .direction = "downup") %>%
  ungroup() %>%
  mutate(
    auc = as.numeric(auc),
    pattern = ifelse(
      is.na(pattern),
      sub(":.*", "", method_pattern),
      pattern
    ),
    methods = ifelse(
      is.na(methods),
      sub(".*:\\s*", "", method_pattern),
      methods
    )
  )
merged_data <- merged_data %>%
  arrange(pattern, methods) %>% 
  mutate(
    method_pattern = factor(method_pattern, levels = unique(method_pattern)),
    pattern_x = method_pattern
  ) %>% 
  arrange(resolution, tech_platform, species, tissue,dt) %>% 
  mutate(dt = factor(dt, levels = unique(dt)))

merged_data <- merged_data %>%
  mutate(pattern = factor(pattern,levels=str_to_sentence(patterns))) %>% 
  arrange(pattern,methods,method_pattern) %>% 
  mutate(method_pattern = factor(method_pattern, levels = unique(method_pattern)))


merged_data <- merged_data %>%
  arrange(resolution, tech_platform, species, tissue, dt) %>% 
  mutate(dt = factor(dt, levels = unique(dt)))

dt_annotations <- merged_data %>%
  distinct(dt, resolution, tech_platform, species, tissue)

pattern_annotations <- merged_data %>%
  distinct(method_pattern, pattern, pattern_x)

palette_mapping <- list(
  resolution = "Set3",       
  tech_platform = "Set2",    
  species = "Paired",        
  tissue = "Dark2",          
  pattern = "Set1"           
)

get_dim_color_palette <- function(categories, palette) {
  n <- length(categories)
  max_col <- brewer.pal.info[palette, "maxcolors"]
  
  if (n > max_col) {
    colors <- colorRampPalette(brewer.pal(max_col, palette))(n)
  } else {
    colors <- brewer.pal(n, palette)
  }
  setNames(colors, categories)
}

res_colors <- get_dim_color_palette(unique(dt_annotations$resolution), palette_mapping$resolution)
tech_colors <- get_dim_color_palette(unique(dt_annotations$tech_platform), palette_mapping$tech_platform)
species_colors <- get_dim_color_palette(unique(dt_annotations$species), palette_mapping$species)
tissue_colors <- get_dim_color_palette(unique(dt_annotations$tissue), palette_mapping$tissue)
pattern_colors <- get_dim_color_palette(unique(pattern_annotations$pattern), palette_mapping$pattern)

annotation_width <- 1
offset_base <- -annotation_width
pattern_anno_height <- 0.6
pattern_offset_base <- -pattern_anno_height
anno_width <- 1
anno_height <- 1
na_df <- merged_data %>% filter(is.na(auc))
na_df <- as.data.frame(na_df)

p1 <- ggplot() +
  geom_tile(
    data = merged_data,
    aes(
      x = method_pattern,
      y = dt,
      fill = auc
    ),
    color = "white",
    width = 1,
    height = 1
  )+
  scale_fill_gradientn(
    colours = colorRampPalette(c("#008bd0", "#eeeeee", "#ffa61d"))(100),
    na.value = "#e0e0e0",   
    values = scales::rescale(c(
      min(merged_data$auc, na.rm = TRUE), 0.7, 0.9 )),
    limits = c(min(merged_data$auc, na.rm = TRUE),0.9),
    oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(0.05, "in"),
      barheight = unit(0.5, "in"),
      label.hjust = 0.5
    ),
    name = "AUC"
  )+
  # resolution
  new_scale_fill() +
  geom_tile(
    data = dt_annotations,
    aes(x = offset_base, y = dt, fill = resolution),
    width = annotation_width, 
    color = NA,  
    linewidth = 0,  
    height = 1,  
    position = "identity"  
  ) +
  scale_fill_manual(values = res_colors, name = "Resolution") +
  guides( 
    fill = guide_legend(
      nrow = 3,           
      byrow = TRUE))+  
  # tech_platform
  new_scale_fill() +
  geom_tile(
    data = dt_annotations,
    aes(x = offset_base - annotation_width, y = dt, fill = tech_platform),
    width = annotation_width, 
    color = NA,  
    linewidth = 0,  
    height = 1,  
    position = "identity"  
  ) +
  scale_fill_manual(values = tech_colors, name = "Tech Platform") +
  guides( 
    fill = guide_legend(nrow = 4, byrow = TRUE))+  
  #  species
  new_scale_fill() +
  geom_tile(
    data = dt_annotations,
    aes(x = offset_base - 2*annotation_width, y = dt, fill = species),
    width = annotation_width, 
    color = NA,  
    linewidth = 0,  
    height = 1,  
    position = "identity" ) +
  scale_fill_manual(values = species_colors, name = "Species") +
  guides( 
    fill = guide_legend(
      nrow = 3,           
      byrow = TRUE))+  
  #  tissue
  new_scale_fill() +
  geom_tile(
    data = dt_annotations,
    aes(x = offset_base - 3*annotation_width, y = dt, fill = tissue),
    width = annotation_width, 
    color = NA,  
    linewidth = 0,  
    height = 1,  
    position = "identity"  
  ) +
  scale_fill_manual(values = tissue_colors, name = "Tissue") +
  guides( fill = guide_legend(ncol = 2,byrow = TRUE))+
  new_scale_fill() +
  geom_tile(
    data = pattern_annotations,
    aes(x = pattern_x, y = pattern_offset_base, fill = pattern),
    height = pattern_anno_height, color = NA
  ) +scale_fill_manual(values = pattern_colors, name = "Pattern") +
  geom_segment(
    data = na_df,
    aes(
      x = as.integer(method_pattern) - 0.4,
      xend = as.integer(method_pattern) + 0.4,
      y = as.integer(dt) - 0.4,
      yend = as.integer(dt) + 0.4
    ),
    colour = "#808080",
    linewidth = 0.25,
    alpha = 0.6
  ) +geom_segment(
    data = na_df,
    aes(
      x = as.integer(method_pattern) - 0.4,
      xend = as.integer(method_pattern) + 0.4,
      y = as.integer(dt) + 0.4,
      yend = as.integer(dt) - 0.4
    ),
    colour = "#808080",
    linewidth = 0.25,
    alpha = 0.6
  )+scale_x_discrete(
    labels = function(x) {
      sub(".*:\\s*", "", x)  
    },
    expand = c(0, 0) 
  ) +
  scale_y_discrete(
    expand = c(0, 0)
  ) +theme_minimal() +
  my_theme+
  theme(
    axis.text.x = element_text(
      angle = 90,   
      hjust = 1,    
      vjust = 0.5
    ),
    axis.title.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.06, "in"),
    legend.key.width = unit(0.06, "in"),
    legend.spacing.x = unit(0.01, "in"),
    legend.spacing.y = unit(0.06, "in"), 
    legend.title = element_text(size = 7,margin = margin(b = 2) ),
    legend.text = element_text(size = 7),
    legend.box = "vertical",
    # legend.box.just = "top",
    legend.margin = margin(1, 1, 0, 0),
    plot.margin = margin(0, 1, 1, 0),  
    panel.grid = element_blank(),     
    panel.border = element_blank()     
  ) +
  labs(y = "Simulated datasets")