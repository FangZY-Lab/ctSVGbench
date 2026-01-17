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
library(cowplot)
source('./my_theme.R')

ncores = 1
files = list.files(here('real','computation'), full.names = TRUE, pattern = 'csv')
files = files[grep("stereoseq_mosta_E16.5_E1S3", files)]
res <- lapply(files, function(f) {
  read.csv(f)
})
res.df = do.call(rbind, res)
res.df$dataset <- recode(res.df$dataset, 'SeqFish+_mouse_ob' = 'SeqFish+_mouse_OB')


res.df <- res.df %>%
  mutate(method = recode(method,
                         "C-SIDE" = "C-SIDE",
                         "CELINA" = "Celina",
                         "STANCE" = "STANCE",
                         "spVC" = "spVC",
                         "ctsvg" = "ctSVG"))

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


# Plot 1: Execution time
p1 <- ggplot(data = res.df) +
  geom_point(aes(x = n_spot, y = time, color = method),size=1) +
  geom_line(aes(x = n_spot, y = time, color = method))+
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
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Run Time (log10 seconds)", fill = "",x="Number of cells") 

tile_height_p2 = 1
CELINA_order_mem <- res.df %>%
  filter(method == "Celina") %>%
  arrange(Peak_RAM_Used_MiB) %>%
  pull(dataset)

res.df$dataset <- factor(res.df$dataset, levels = CELINA_order_mem)
p2 <- ggplot(res.df) +
  geom_point(aes(x = n_spot, y = Peak_RAM_Used_MiB/1024, color = method), size = 1) +
  geom_line(aes(x = n_spot, y = Peak_RAM_Used_MiB/1024, color = method)) +
  scale_color_manual(values = method_colors) +
  scale_y_log10(
    breaks = c(0.1, 1, 10, 20,30),
    labels = scales::label_number()
  )  +
  theme_minimal() +
  my_theme +
  theme(
    legend.position = 'left',
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Memory (GiB)", color = "Methods",x="Number of cells") 
p2  

plot_grid(p1,  p2, rel_widths = c(1,1.02), ncol = 2, labels = c("A", "B"))
ggsave("./Fig/Fig5.pdf", width = 6.69, height = 3)


# Summary statistics
summary_df_sca <- res.df[(res.df$n_spot>1500&res.df$n_spot<2500),] %>%
  group_by(method) %>%
  summarise(
    CPU_time = median(time / 3600, na.rm = TRUE),
    RAM_GB = median(Peak_RAM_Used_MiB / 1024, na.rm = TRUE),
  ) %>% column_to_rownames('method')

summary_df <- read.csv("metrics_summary.csv", row.names = 1)
summary_df$CPU_time_2000 <- summary_df_sca[rownames(summary_df),"CPU_time"]
summary_df$RAM_GB_2000 <- summary_df_sca[rownames(summary_df),"RAM_GB"]
summary_df_sca <- res.df[(res.df$n_spot>4500&res.df$n_spot<5500),] %>%
  group_by(method) %>%
  summarise(
    CPU_time = median(time / 3600, na.rm = TRUE),
    RAM_GB = median(Peak_RAM_Used_MiB / 1024, na.rm = TRUE),
  ) %>% column_to_rownames('method')

summary_df$CPU_time_5000 <- summary_df_sca[rownames(summary_df),"CPU_time"]
summary_df$RAM_GB_5000 <- summary_df_sca[rownames(summary_df),"RAM_GB"]

write.csv(summary_df, "metrics_summary.csv", row.names = TRUE)
summary_df$CPU_time<- NULL
summary_df$RAM_GB <-NULL
