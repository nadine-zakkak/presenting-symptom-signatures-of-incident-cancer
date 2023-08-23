## --------------------------------
## Script name: 3_figures_tables.R
## 
## Previous Script: 2_analysis.R
##
## Purpose of Script: Produce all required tables and figures for the paper
##
## Author: Nadine Zakkak
##
## ---------------------------------

# load libraries needed for all outputs in the beginning -----
require('pacman')
pacman::p_load(plyr, dplyr, tidyr, stringr,
               gt, ggplot2, ggsci, ggnewscale)

# output directory
output_dir <- "./results/final-tables-figs/"

# load functions and global var
source("./scripts/00_global_var.R")
source("./scripts/00_functions.R")

# Table 1 - Sample characteristics by cancer site -----
table1 <- readRDS("./data/table1_info.rds")

table1 |>
  mutate(cancer_group = factor(cancer_group, levels = c("Total", cancergrp_order)),
         cancer_site_desc = factor(cancer_site_desc, levels = c("Total", cancer_order))) |>
  arrange(cancer_group, cancer_site_desc) |>
  mutate(across(everything(), ~replace_na(as.character(.), "-"))) |>
  select(cancer_group, cancer_site_desc, 
         `N_1`, `summary_1_<60`, `summary_1_60-74`, `summary_1_75+`,
         `N_2`, `summary_2_<60`, `summary_2_60-74`, `summary_2_75+`) |>
  rename(`Cancer Site` = cancer_site_desc) |>
  group_by(cancer_group) |>
  gt() |>
  tab_spanner(label = "Men",
              columns = c(3,4,5,6)) |>
  tab_spanner(label = "Women",
              columns = c(7,8,9,10)) |>
  cols_label(
    N_1 = "N",
    `summary_1_<60` = "<60",
    `summary_1_60-74` =   "60-74",
    `summary_1_75+` = "75+",
    N_2 = "N",
    `summary_2_<60` = "<60",
    `summary_2_60-74` =   "60-74",
    `summary_2_75+` = "75+"
  ) |>
  gtsave(paste0(output_dir, "table1.rtf"))

rm(table1)

# Table 2 - Number of symptoms per cancer site that occurred in >1% and >50% of patients ----
numb_sx_perc <- readRDS("./data/numb_sx_perc.rds")

numb_sx_perc |>
  select(-c(`>5%`, `>10%`, `>20%`)) |>
  mutate(cancer_site_desc = factor(cancer_site_desc, levels = cancer_order)) |>
  arrange(cancer_site_desc) |>
  rename("Cancer Site" = cancer_site_desc) |>
  group_by(cancer_group) |>
  gt() |>
  tab_spanner(label = "Number of symptoms in", columns = matches(">.*")) |>
  gtsave(paste0(output_dir, "table2.rtf"))

rm(numb_sx_perc)

# Figure 2 - Heatmap at group level (proportion of symptom group in each cancer group) ----
htmap_grp_df <- readRDS("./data/htmap_grp_df.rds")
htmap_grp_txt_col <- readRDS("./data/htmap_grp_txt_col.rds")

htmap_col <- c('low' = "#FFFEF8", 'high' = '#705BD0')

htmap_grp <- 
  htmap_grp_df |>
  ggplot(aes(x = cancer_group, 
             y = symptom_group, 
             fill = prop)) +
  geom_tile(data = htmap_grp_df,
            color = "white") +
  scale_fill_gradient(name = "Column Percentage",
                      low = htmap_col['low'], 
                      high = htmap_col['high'], 
                      limits = c(0, 99)) +
  geom_text(aes(label = ifelse(prop == 0, "-", prop)), size = 6, color = htmap_grp_txt_col) +
  scale_x_discrete(position = "top") +
  ylim(rev(sxgrp_order[sxgrp_order != "summary"])) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 1, size = 15, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) +
  labs(x = "Cancer Groups", y = "Symptom Groups")

ggsave(plot = htmap_grp, filename  = paste0(output_dir, "figure2.png"),
       units = "cm", width = 35, height = 30)

rm(list = ls(pattern = "^htmap"))

# Figure 3 - Heatmap at individual level (proportion of symptom in each cancer site) ----
htmap_df <- readRDS("./data/htmap_df.rds")
htmap_cancersite_col <- readRDS("./data/htmap_cancersite_col.rds")
htmap_sx_col <- readRDS("./data/htmap_sx_col.rds")

htmap_col <- c('low' = "#FFFEF8", 'high' = '#705BD0')

htmap_df$cancer_group <- factor(htmap_df$cancer_group, levels = cancergrp_order)
htmap_df$cancer_site_desc <- factor(htmap_df$cancer_site_desc, levels = cancer_order)

htmap_df$symptom_group <- factor(htmap_df$symptom_group, levels = sxgrp_order)
htmap_df$symptom_fact <- factor(htmap_df$symptom_recod, 
                                levels = unique(htmap_df$symptom_recod[order(htmap_df$symptom_group)]))

htmap_vline <- # helper df for legend groups
  htmap_df |> 
  distinct(cancer_group, cancer_site_desc) |> 
  group_by(cancer_group) |> # 1 row per cancer group
  summarise(n = n()) |>
  mutate(xintercept = cumsum(n) + .5) # coord to be used to draw vertical lines

htmap_hline <-
  htmap_df |>
  distinct(symptom_group, symptom_fact) |>
  group_by(symptom_group) |> # 1 row per cancer group
  summarise(n = n()) |>
  mutate(y = cumsum(n) + .5,
         yintercept = 85 - y) # coord to be used to draw horizontal lines

htmap_txt_col <- ifelse(round(htmap_df$prop, 0) == 0, "#d3d3d3", "black")

htmap <- htmap_df |>
  ggplot(aes(x = cancer_site_desc, 
             y = symptom_fact, 
             fill = prop)) +
  geom_tile(data = filter(htmap_df, symptom_fact != 'Mean no. symptoms'), 
            color = "white") +
  scale_fill_gradient(name = "Column Percentage",
                      low = htmap_col['low'], 
                      high = htmap_col['high'], 
                      limits = c(0, 99)) +
  new_scale_fill() +
  geom_tile(data = filter(htmap_df, symptom_fact == 'Mean no. symptoms'), 
            mapping = aes(fill = symptom_fact == 'Mean no. symptoms'),
            color = "white") +
  scale_fill_manual(values = "grey", 
                    guide = "none",
                    labels = "") +
  geom_text(aes(label = ifelse(prop == 0, "-", prop)), size = 18/.pt, color = htmap_txt_col) +
  geom_vline(xintercept = htmap_vline$xintercept, linewidth = .1, color = "#aeaeae") +
  geom_hline(yintercept = htmap_hline$yintercept, linewidth = .1, color = '#aeaeae') +
  scale_x_discrete(position = "top") +
  ylim(rev(levels(htmap_df$symptom_fact))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 1, size = 16, color = htmap_cancersite_col),
        axis.text.y = element_text(size = 16, color = htmap_sx_col),
        axis.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom") +
  labs(x = "Cancer Sites", y = "Symptoms")

ggsave(plot = htmap, filename  = paste0(output_dir, "figure3.png"),
       units = "cm", width = 60, height = 60)

rm(list = ls(pattern = "^htmap"))

# Figure 4 - Barplot at group level -----
colours <- c("#EE5250", "#AC1357", "#B967C7", "#FF6E00",
             "#C5CAE9", "#0C46A0", "#02A9F3", "#1A5E1F",
             "#DCECC7", "#004C3F", "#FFEBED", "#817717",
             "#A0877F", "#C4B3AB", "#FFF2DF", "#202020")

barplt_grp_df <- readRDS("./data/barplt_grp_df.rds")

barplt_group <- 
  barplt_grp_df |>
  ggplot(aes(x = symptom_group, y = prop, fill = cancer_group)) +
  geom_bar(stat = "identity", aes(fill = cancer_group)) +
  scale_fill_manual(values = colours, name = "Cancer Groups") +
  facet_grid(~sex, scales = "free") +
  coord_flip() +
  labs(x = "", 
       y = "") + 
  theme_minimal() +
  theme(strip.text.x = element_text(size = 13, color = "black"),
        text = element_text(size = 13, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.key.size = unit(1.5, "lines"))

ggsave(plot = barplt_group, filename = paste0(output_dir, "figure4.png"),
       units = "cm", width = 35, height = 30)

rm(list = ls(pattern = "^barplt"))


# Figure 5 - Barplot at individual level -----
sxgrp_fig_name <- read.delim("./lookup_files/sxgrp_figures_rename.txt")
barplt_df <- readRDS("./data/barplt_df.rds")
barplt_legend_grp <- readRDS("./data/barplt_legend_grp.rds")
barplt_labels <- readRDS("./data/barplt_labels.rds")
barplt_colors <- readRDS("./data/barplt_colors.rds")

barplt_df$cancer_group <- factor(barplt_df$cancer_group, levels = cancergrp_order)
barplt_df$cancer_site_desc <- factor(barplt_df$cancer_site_desc, levels = cancer_order)

barplt_df <- barplt_df |> arrange(desc(symptom_recod))
barplt_df$symptom_group <- factor(barplt_df$symptom_group, levels = unique(mapvalues(sxgrp_order,
                                                                                     from = sxgrp_fig_name$Old,
                                                                                     to = sxgrp_fig_name$New)))
barplt_df$symptom_recod <- factor(barplt_df$symptom_recod, 
                                  levels = unique(barplt_df$symptom_recod[order(barplt_df$symptom_group)]))
add_legend <- function(i){ # helper function to create legends
  do.call(scale_fill_manual, 
          list(aesthetics = "fill", 
               values = barplt_colors, 
               labels = barplt_labels[from[i]:to[i]], 
               breaks = names(barplt_colors)[from[i]:to[i]], 
               name = cancer_group[i], 
               guide = guide_legend(order = i, ncol = 1)))
}

attach(barplt_legend_grp)
{ # plot ----
  p <- barplt_df |> 
    ggplot(aes(x = symptom_recod, y = prop, fill = cancer_site_desc)) +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(1) + 
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(2) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(3) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(4) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(5) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(6) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(7) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(8) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(9) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(10) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(11) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(12) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(13) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(14) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(15) +
    new_scale_fill() +
    geom_bar(stat = "identity", aes(fill = cancer_site_desc)) +
    add_legend(16) +
    facet_grid(symptom_group~sex, scales = "free", space = "free") +
    coord_flip() +
    labs(x = "", 
         y = "") + 
    theme_minimal() +
    theme(strip.text.x = element_text(size = 13, colour = "black"),
          strip.text.y = element_text(size = 10.5, colour = "black"),
          text = element_text(size = 11, colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 13, colour = "black"),
          legend.text = element_text(size = 13, colour = "black"),
          legend.title = element_text(size = 13, colour = "black"),
          panel.spacing.y = unit(1, "lines")
    )
}
detach(barplt_legend_grp)

ggsave(plot = p, filename = paste0(output_dir, "figure5.png"),
       units = "cm", width = 45, height = 50)

rm(sxgrp_fig_name, p, list = ls(pattern = "^barplt"))
