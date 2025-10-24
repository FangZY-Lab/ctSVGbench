library(devtools)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(ggrepel)
library(here)
library(ggnewscale)
library(RColorBrewer)
library(tidyverse)
source('F:/ctSVGbench/my_theme.R')

ncores = 1
files = list.files(here('real','computation'), full.names = TRUE, pattern = 'csv')
files = files[-grep("sim", files)]
files = files[-grep("CTSV", files)]
res <- lapply(files, function(f) {
  read.csv(f)
})
res.df = do.call(rbind, res)
res.df$dataset <- recode(res.df$dataset, 'SeqFish+_mouse_ob' = 'SeqFish+_mouse_OB')

dts <- c(
  "MERFISH_hypothalamus",
  "SeqFish+_mouse_OB",
  "ST_PDAC",
  "StereoSeq_CBMSTA_Marmoset",
  "Visium_liver",
  "Visium_skin",
  "Visium_spleen",
  "Visium_bladder",
  "Visium_tail",
  "StereoSeq_mouseOB",
  "Visium_mousebrain",
  "ST_developmental_heart",
  "Visium_intestine",
  "Visium_pancreas",
  "StereoSeq_MDESTA",
  "Slide-seqV2_melanoma",
  "Visium_lymph_node",
  "Slide-seq_tumor",
  "Slide-seqV2_hippocampus",
  "StereoSeq_CBMSTA_Macaque",
  "Slide-seqV2_mouseOB",
  "Visium_melanoma"
)  

methods <- c("C-SIDE", "spVC", "CELINA", "STANCE") 
files <- character()  
for (mt in methods){
  for (dt in dts) {
    file <- sprintf("%s-%s.csv", dt, mt)
    files <- c(files, file)  
  }
}

res <- lapply(files, function(f){
  read.csv(here('real','computation', f))
})
res.df = do.call(rbind, res)

res.df <- res.df %>%
  mutate(method = recode(method,
                         "C-SIDE" = "C-SIDE",
                         "CELINA" = "Celina",
                         "STANCE" = "STANCE",
                         "spVC" = "spVC"))

# Create label information
res.df <- res.df %>%
  mutate(label_info = paste0(n_spot, " spots\n", n_gene, " genes"))
res.df$dataset <- recode(res.df$dataset, 'SeqFish+_mouse_ob' = 'SeqFish+_mouse_OB')

# Convert seconds to formatted time (days/hours/minutes)
sec_to_time <- function(x) {
  d <- floor(x / 86400)
  h <- floor((x %% 86400) / 3600)
  m <- floor((x %% 3600) / 60)
  s <- round(x %% 60)
  
  label <- ifelse(d > 0, sprintf("%dd %02dh", d, h),
                  ifelse(h > 0, sprintf("%dh %02dm", h, m),
                         sprintf("%dm %02ds", m, s)))
  return(label)
}

log_breaks <- c(1, 10, 100, 1000, 10000, 100000)
library(ggnewscale)
CELINA_order_time <- res.df %>%
  filter(method == "Celina") %>%
  arrange(time) %>%
  pull(dataset)

res.df$dataset <- factor(res.df$dataset, levels = CELINA_order_time)

res.df <- res.df %>%
  mutate(label_info = paste0(n_spot, " spots\n", n_gene, " genes"))

tile_height_p1 <- 0.02 * (max(res.df$time) - min(res.df$time))

# Plot 1: Execution time
p1 <- ggplot(data = res.df) +
  geom_point(aes(x = dataset, y = time, color = method, size = log10(time))) +
  scale_size_continuous(range = c(0.2, 4)) +
  scale_color_manual(values = method_colors) +
  theme_minimal() +
  my_theme +
  scale_y_log10(
    breaks = log_breaks,
    labels = log_breaks,
    sec.axis = sec_axis(
      ~ .,
      name = "Run Time (Days/Hours/Minutes)",
      breaks = log_breaks,
      labels = sec_to_time
    )
  ) +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Run Time (log10 seconds)", fill = "") +
  new_scale_fill() +
  geom_tile(aes(x = dataset, y = 0.7, fill = n_gene),
            height = 0.1 * (log10(max(res.df$time)) - log10(min(res.df$time))),
            width = 1) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "n_gene") +
  new_scale_fill() +
  geom_tile(aes(x = dataset, y = 0.9, fill = n_spot),
            height = 0.22, width = 1) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu")), name = "n_spot")

# Plot 2: Memory usage
tile_height_p2 = 1
CELINA_order_mem <- res.df %>%
  filter(method == "Celina") %>%
  arrange(Peak_RAM_Used_MiB) %>%
  pull(dataset)

res.df$dataset <- factor(res.df$dataset, levels = CELINA_order_mem)

p2 <- ggplot(res.df) +
  geom_point(aes(x = dataset, y = Peak_RAM_Used_MiB/1024, size = Peak_RAM_Used_MiB, color = method)) +
  scale_size_continuous(range = c(0.2, 4)) +
  scale_color_manual(values = method_colors) +
  theme_minimal() +
  my_theme +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    plot.margin = margin(10, 49, 10, 22)
  ) +
  ylim(-2.5, (max(res.df$Peak_RAM_Used_MiB/1024) + 1)) +
  labs(y = "Memory (GiB)", fill = "") +
  new_scale_fill() +
  geom_tile(aes(x = dataset, y = -1.0, fill = n_spot), height = tile_height_p2, width = 1) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu")), name = "n_spot") +
  new_scale_fill() +
  geom_tile(aes(x = dataset, y = -2.0, fill = n_gene), height = tile_height_p2, width = 1) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "n_gene")

p2

# Plot 3: Legend extraction
p3 <- ggplot(res.df, aes(x = dataset, y = time, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  scale_fill_manual(values = method_colors) +
  theme_void() +
  my_theme +
  theme(
    plot.margin = margin(0, 0, 0, 5),
    legend.position = 'bottom',
    legend.key.size = unit(0.1, "in"), 
    legend.key.width = unit(0.1, "in")
  ) +
  labs(fill = '') 
ggsave('p3.pdf')

l3 <- ggpubr::get_legend(p3)

# Plot 4: Legend tiles
p4 <- res.df %>% ggplot() +
  new_scale_fill() +
  geom_tile(aes(x = dataset, y = 0.9, fill = n_gene), height = 0.05, width = 1) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"),
                       name = "n_gene",
                       breaks = c(min(res.df$n_gene), median(res.df$n_gene), max(res.df$n_gene)),
                       labels = scales::comma,
                       guide = guide_colorbar(
                         title.position = "top",
                         barwidth = unit(2.5, "cm"),
                         barheight = unit(0.2, "cm"),
                         label.hjust = 0.5
                       )) +
  new_scale_fill() +
  geom_tile(aes(x = dataset, y = 0.7, fill = n_spot), height = 0.05, width = 1) +
  scale_fill_gradientn(
    colors = rev(brewer.pal(11, "RdYlBu")), name = "n_spot",
    breaks = c(min(res.df$n_spot), median(res.df$n_spot), max(res.df$n_spot)),
    labels = scales::comma,
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(2.5, "cm"),
      barheight = unit(0.2, "cm"),
      label.hjust = 0.5
    )
  ) +
  theme_void() +
  my_theme +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(1, "cm"),
    legend.key.size = unit(0.5, "in"), 
    legend.key.width = unit(0.2, "in"), 
    plot.margin = margin(0, 120, 0, 0)
  )

p4
ggsave("p4.pdf", p4, width = 6, height = 4)

l4 <- ggpubr::get_legend(p4)

library(cowplot)
l34 <- plot_grid(l3, l4, rel_widths = c(1, 1.2), ncol = 2)
plot_grid(p1, l34, p2, rel_heights = c(1, 0.2, 1), nrow = 3, labels = c("A", "", "B"))
ggsave("./Fig/Fig5-2.pdf", width = 6.69, height = 8)

# Summary statistics
summary_df_sca <- res.df %>%
  group_by(method) %>%
  summarise(
    CPU_time = median(time / 3600, na.rm = TRUE),
    RAM_GB = median(Peak_RAM_Used_MiB / 1024, na.rm = TRUE),
  ) %>% column_to_rownames('method')

summary_df <- read.csv("metrics_summary.csv", row.names = 1)
summary_df$CPU_time <- summary_df_sca[rownames(summary_df),"CPU_time"]
summary_df$RAM_GB <- summary_df_sca[rownames(summary_df),"RAM_GB"]
write.csv(summary_df, "metrics_summary.csv", row.names = TRUE)
