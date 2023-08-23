## --------------------------------
## Script name: 4_supplementary_material.R
## 
## Previous Script: 3_figures_tables.R
##
## Purpose of Script: Produce all required supplementary material for the paper
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

# Supplementary table 1 - ICD-10, cancer site and cancer group classification -----
final_cancer_class <- readRDS("./data/final_cancer_class.rds")

final_cancer_class |>
  select(icd10, cancer_group, cancer_site_desc) |>
  rename("ICD10 O2"=icd10,
         "Cancer Group" = cancer_group,
         "Cancer Site" = cancer_site_desc) |>
  gt() |>
  gtsave(paste0(output_dir, "supp1.rtf"))

# Supplementary table 2 - Symptoms classification -----
final_symptom_class <- readRDS("./data/final_symptom_class.RDS")

final_symptom_class |>
  gt() |>
  gtsave(paste0(output_dir, "supp2.rtf"))

# Supplementary table 3 - Proportion of symptom group by cancer group  -------
sx_site_grp_table <- readRDS("./data/sx_site_grp_table.rds")

sx_site_grp_table <- 
  sx_site_grp_table |>
  arrange(cancer_group, symptom_group) |>
  mutate(`Cancer Group` = str_glue("{cancer_group} (N = {N_site})")) |>
  select(-c(N_site, cancer_group)) |>
  rename("% (95% CI)" = summary)

sx_site_grp_table |>
  group_by(`Cancer Group`) |>
  gt() |>
  gtsave(paste0(output_dir, "supp3.rtf"))

write.csv(sx_site_grp_table |>
            select(`Cancer Group`, symptom_group, everything())
          ,file = paste0(output_dir, "supp3.csv"), row.names = F)

rm(sx_site_grp_table)

# Supplementary table 4 - Proportion of symptom by cancer site ------
sx_site_table <- readRDS("./data/sx_site_table.rds")

sx_site_table <- 
  sx_site_table |>
  arrange(symptom) |>
  mutate(symptom = factor(symptom, levels = unique(symptom[order(symptom_group)])),
         cancer_group = factor(cancer_group, levels = cancergrp_order),
         cancer_site_desc = factor(cancer_site_desc, levels = cancer_order)) |>
  arrange(cancer_site_desc, symptom) |>
  select(-cancer_site_desc) |>
  rename("% (95% CI)" = summary)

sx_site_table |>
  select(symptom_group, everything()) |>
  group_by(cancer_group, `Cancer Site`) |>
  gt()  |>
  gtsave(paste0(output_dir, "supp4.rtf"))

write.csv(sx_site_table |>
            select(symptom_group, symptom, `Cancer Site`, everything())
          ,file = paste0(output_dir, "supp4.rtf"), row.names = F)

rm(sx_site_table)

# Supplemenatry table 5 - Number of symptoms per cancer site that occurred in >1%, >5%, >10%, >20% and >50% of patients ------
numb_sx_perc <- readRDS("./data/numb_sx_perc.rds")

numb_sx_perc |>
  mutate(cancer_site_desc = factor(cancer_site_desc, levels = cancer_order)) |>
  arrange(cancer_site_desc) |>
  rename("Cancer Site" = cancer_site_desc) |>
  group_by(cancer_group) |>
  gt() |>
  tab_spanner(label = "Number of symptoms in", columns = matches(">.*")) |>
  gtsave(paste0(output_dir, "supp5.rtf"))

rm(numb_sx_perc)

# Supplementary table 6 - Proportions and CI for cancer group by sx groups (stratified by sex) -----
site_sx_grp_table <- readRDS("./data/site_sx_grp_table.rds")

site_sx_grp_table <- 
  site_sx_grp_table |>
  arrange(symptom_group, cancer_group) |>
  rename("% (95% CI)" = summary,
         `Cancer Group` = cancer_group) |>
  pivot_wider(id_cols = c(symptom_group, `Cancer Group`),
              names_from = sex,
              values_from = c(N_sx, `% (95% CI)`)) |>
  group_by(symptom_group) |>
  fill(N_sx_Women, .direction = "down") |>
  fill(N_sx_Men, .direction = "down") |>
  ungroup() |>
  mutate(across(everything(), ~replace_na(as.character(.), "-")), 
         Symptom = str_glue("{symptom_group} (Men={N_sx_Men}, Women={N_sx_Women})"))

site_sx_grp_table |>
  select(Symptom, `Cancer Group`, matches("\\%.*Men$"), matches("\\%.*Women$")) |>
  group_by(Symptom) |>
  gt() |>
  tab_spanner("Men", columns = 3) |>
  tab_spanner("Women", columns = 4) |>
  gtsave(paste0(output_dir, "supp6.rtf"))

write.csv(site_sx_grp_table |>
            select(-Symptom)
          ,file = paste0(output_dir, "supp6.csv"), row.names = F)

rm(site_sx_grp_table)

# Supplementary table 7 - Proportions and CI for site by sx (stratified by sex) -----
site_sx_table <- readRDS("./data/site_sx_table.rds")

site_sx_table <- 
  site_sx_table |>
  arrange(symptom, cancer_site_desc) |>
  rename("% (95% CI)" = summary,
         `Cancer Site` = cancer_site_desc,
         `Cancer Group` = cancer_group,
         `Symptom Group` = symptom_group) |>
  pivot_wider(id_cols = c(`Symptom Group`, symptom, `Cancer Group`, `Cancer Site`),
              names_from = sex,
              values_from = c(N, `% (95% CI)`)) |>
  group_by(`Symptom Group`) |>
  fill(N_Women, .direction = "down") |>
  fill(N_Men, .direction = "down") |>
  ungroup() |>
  mutate(across(everything(), ~replace_na(as.character(.), "-")), 
         Symptom = str_glue("{`Symptom Group`} - {symptom} (Men={N_Men}, Women={N_Women})"))

site_sx_table |>
  select(Symptom, `Cancer Group`, `Cancer Site`, matches("\\%")) |>
  group_by(Symptom) |>
  gt() |>
  tab_spanner("% (95% CI)", columns = c(4, 5)) |>
  gtsave(paste0(output_dir, "supp7.rtf"))

write.csv(site_sx_table |>
            select(`Symptom Group`, Symptom, `Cancer Group`, `Cancer Site`,
                   matches(".*_Men"), matches(".*_Women"))
          ,file = paste0(output_dir, "supp7.csv"), row.names = F)

rm(site_sx_table)
