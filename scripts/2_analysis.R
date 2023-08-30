## --------------------------------
## Script name: 2_analysis.R
## 
## Previous Script: 1_preprocessing.R
##
## Purpose of Script: Analysis of NCDA data for symptom signature study
##
## Author: Nadine Zakkak
##
## ---------------------------------

# Loading libraries ----------------------
require('pacman')
pacman::p_load(readr, plyr, dplyr, ggplot2, tidyr, fastDummies, 
               gt, stringr, forcats, ggsci, ggnewscale)

# Load data into R ------
# data <- readRDS("/path/to/data/single_tumour_wide.RDS")
# data_long <- readRDS("/path/to/data/single_tumour_long.RDS")
# cohort <- readRDS("/path/to/data/cohort.RDS")
cancer_groups <- read.delim("./lookup_files/cancer_site_groups.txt") 
symptom_groups <- read.delim("./lookup_files/symptom_groups.txt", sep=" ")


# Location to save files ----
output_dir <- "./results/"

# Load functions and global variales
source("./scripts/00_global_var.R")
source("./scripts/00_functions.R")

# ICD - Cancer site from the data ---- 
cancer_codelist_data <- data |> distinct(site_icd10_o2, cancer_site_desc, cancer_group)
cancer_codelist_data <-
  cancer_codelist_data |>
  mutate(icd10_3dig = substr(site_icd10_o2, 1, 3)) |>
  arrange(site_icd10_o2) |>
  group_by(cancer_site_desc, cancer_group) |>
  summarise(icd10_4dig = toString(site_icd10_o2),
            icd10_3dig = toString(icd10_3dig)) |>
  rowwise() |>
  mutate(icd10_3dig = toString(unique(unlist(strsplit(icd10_3dig, ", ")))))

# Cohort - general description ----
cohort |>
  group_by(sex) |>
  summarise(n = n()) |>
  mutate(prop = 100*calculate_prop(n, sum(n)),
         covariate = "sex",
         sex = ifelse(sex==1, "Men", "Women")) |>
  rename(value = sex) |>
  bind_rows(
    cohort |>
      group_by(fiveyearageband) |>
      summarise(n = n()) |>
      mutate(prop = 100*calculate_prop(n, sum(n)),
             covariate = "five year age band") |>
      rename(value = fiveyearageband)) |>
  bind_rows(
    cohort |>
      group_by(tenyearageband) |>
      summarise(n = n()) |>
      mutate(prop = 100*calculate_prop(n, sum(n)),
             covariate = "ten year age band") |>
      rename(value = tenyearageband)) |>
  bind_rows(cohort |>
              mutate(var = case_when(
                fiveyearageband %in% c("25-29", "30-34", "35-39",
                                       "40-44", "45-49", "50-54", "55-59") ~ "<60",
                fiveyearageband %in% c("60-64", "65-69", "70-74") ~ "60-74",
                fiveyearageband %in% c("75-79", "80-84", "85-89", "90+") ~ "75+")) |>
              group_by(var) |>
              summarise(n = n()) |>
              mutate(prop = 100*calculate_prop(n, sum(n)),
                     covariate = "3 age groups") |>
              rename(value = var)) |>
  group_by(covariate) |>
  mutate(prop = round(prop, 2)) |>
  gt() |>
  gtsave(paste0(output_dir, "cohort_characterisation.rtf"))


# Useful dataframes (site/symptom frequencies) -----
{
  ## Cancer site frequency ----
  site_freq <-
    data |>
    group_by(cancer_site_desc, cancer_group) |>
    summarise(N = n()) |>
    arrange(desc(N)) |>
    ungroup()
  
  ## Symptom frequency ----
  sx_freq <-
    data_long |>
    group_by(symptom, symptom_group) |>
    summarise(N_sx = sum(present)) |>
    arrange(desc(N_sx)) |>
    ungroup()
  
  ## % of symptom in X cancer site ----
  sx_site <-
    data_long |>
    group_by(cancer_site_desc, symptom) |>
    summarise(N_site = n(),
              n = sum(present),
              prop = mean(present)*100) |>
    ungroup() |>
    inner_join(symptom_groups) |>
    inner_join(cancer_codelist_data |> select(-icd10_3dig, -icd10_4dig))
  
  ## % of X cancer site in symptom Y ----
  site_sx <- data_long |>
    group_by(symptom, cancer_site_desc) |>
    summarise(n = sum(present)) |>
    mutate(N_sx = sum(n),
           prop = 100*calculate_prop(n, N_sx)) |>
    ungroup() |>
    inner_join(symptom_groups) |>
    inner_join(cancer_codelist_data |> select(-icd10_3dig, -icd10_4dig))
  
  ## % of symptom group in X cancer group ----
  sx_site_group <-
    data_long |>
    group_by(patient_pseudoid, tumour_pseudoid, cancer_group, symptom_group) |>
    summarise(count = max(present)) |>
    ungroup() |>
    group_by(cancer_group, symptom_group) |>
    summarise(N_site = n(), # number of patients in each cancer group
              n = sum(count)) |> # number of patients in cancer group x symptom group 
    mutate(prop = 100*calculate_prop(n, N_site)) |>
    ungroup() |>
    mutate(cancer_group = factor(cancer_group, levels = cancergrp_order),
           symptom_group = factor(symptom_group, levels = sxgrp_order)) |>
    arrange(cancer_group, symptom_group)
  
  ## % of X cancer group in symptom group (stratified by sex) ----
  site_sx_group <- 
    data_long |>
    group_by(sex, symptom_group, cancer_group) |>
    summarise(n = sum(present)) |>
    mutate(N_sx = sum(n),
           prop = 100*calculate_prop(n, N_sx)) |>
    ungroup() |>
    mutate(cancer_group = factor(cancer_group, levels = cancergrp_order),
           symptom_group = factor(symptom_group, levels = sxgrp_order)) |>
    arrange(sex, symptom_group, cancer_group)
  
  ## Symptom frequency by sex ----
  sx_sex <-
    data_long |> 
    group_by(sex, symptom) |> 
    mutate(count = 1) |> 
    summarise(N_sex = sum(count),
              n = sum(present),
              prop = mean(present)*100, 1,
              summary = str_glue("{n} ({format_numb(prop)}%)")) |> 
    ungroup() |>
    inner_join(symptom_groups) |>
    mutate(symptom_group = factor(symptom_group, levels=sxgrp_order)) |>
    arrange(symptom_group, symptom) 
  
  ## Site frequency by sex ----
  site_sex <- data |> 
    group_by(sex, cancer_site_desc) |> 
    summarise(n = n()) |> 
    mutate(N = sum(n), 
           prop = 100*calculate_prop(n, N),
           summary = str_glue("{n} ({format_numb(prop)}%)")) |> 
    ungroup() |>
    inner_join(cancer_codelist_data) |>
    mutate(cancer_site_desc = factor(cancer_site_desc, levels = cancer_order)) |>
    arrange(cancer_site_desc)
  
  ## tidying all in a list ----
  freq_tables <- list('site_freq' = site_freq, 'sx_freq' = sx_freq, 
                      'sx_site' = sx_site, 'site_sx' = site_sx,
                      'sx_site_group' = sx_site_group,
                      'site_sx_group' = site_sx_group,
                      'sx_sex' = sx_sex, 
                      'site_sex' = site_sex)
  
  rm(site_freq, sx_freq, sx_site, site_sx, sx_site_group, site_sx_group, sx_sex, site_sex)
}


# Tables - Symptom & Site frequency by sex ----
# ## Site frequency by sex ----
site_sex_wide <- 
  freq_tables$site_sex |>
  pivot_wider(id_cols = c(cancer_site_desc, cancer_group, icd10_4dig, icd10_3dig),
              names_from = c(sex),
              values_from = c(n, prop, summary))
site_sex_wide <-
  site_sex_wide |>
  mutate(total_n = rowSums(site_sex_wide[, c('n_1', 'n_2')], na.rm = T)) |>
  add_row(cancer_site_desc = 'Total', 
          summary_1 = str_glue("{colSums(site_sex_wide[, 'n_1'], na.rm = T)}"),
          summary_2 = str_glue("{colSums(site_sex_wide[, 'n_2'], na.rm = T)}"))

site_sex_wide |>
  rename(`Cancer Site` = cancer_site_desc,
         `ICD-10 (4 dig)`  = icd10_4dig,
         `ICD-10 (3 dig)` = icd10_3dig,
         Men = "summary_1",
         Women = "summary_2") |>
  group_by(cancer_group) |>
  select(-c(n_1, n_2, prop_1,prop_2, total_n)) |> 
  gt() |>
  cols_move_to_end(c(`ICD-10 (4 dig)`, `ICD-10 (3 dig)`)) |>
  gtsave(paste0(output_dir, "site_frequency_bysex.rtf"))

write.csv(site_sex_wide, paste0(output_dir, "site_frequency_bysex.csv"))

rm(site_sex_wide)

## Symptom frequency by sex ----
sx_sex_wide <- 
  freq_tables$sx_sex |>
  mutate(n = ifelse(n<20, "N<20", n),
         summary = ifelse(n == "N<20", "-", summary)) |>
  mutate(prop = format_numb(prop)) |>
  pivot_wider(id_cols = c(symptom, symptom_group),
              names_from = c(sex),
              values_from = c(n, prop, summary))

sx_sex_wide |>
  rename(Men = "summary_1",
         Women = "summary_2") |>
  group_by(symptom_group) |>
  select(-c(n_1, n_2, prop_1,prop_2)) |>
  gt() |>
  gtsave(paste0(output_dir, "symptom_frequency_bysex.rtf"))

write.csv(sx_sex_wide, paste0(output_dir, "symptom_frequency_bysex.csv"))

rm(sx_sex_wide)


# Tables - Missing Symptoms ----
## Missing Symptoms by site ----
# 'row percentages'
freq_tables$sx_site |>
  filter(symptom_group == "None recorded") |>
  mutate(prop = round(prop, 1)) |>
  mutate(n = ifelse(n<20, "n<20", n),
         prop = ifelse(n == "n<20", "-", prop)) |>
  arrange(desc(N_site)) |>
  pivot_wider(id_cols = c(cancer_site_desc, N_site), 
              names_from = symptom, 
              values_from = c('n', 'prop')) |>
  mutate(cancer_site_desc = str_glue("{cancer_site_desc} (N = {N_site})")) |>
  select(-N_site) |>
  gt() |>
  tab_spanner(columns = c("n_N/A", "prop_N/A"), label = "N/A") |>
  tab_spanner(columns = c("n_N/K", "prop_N/K"), label = "N/K") |>
  cols_label("n_N/A" = "n",
             "prop_N/A" = "%",
             "n_N/K" = "n",
             "prop_N/K" = "%",
             "cancer_site_desc" = "Cancer Site") |>
  gtsave(paste0(output_dir, "missing_symptom_by_site.rtf"))

## Missing Symptom overall ----
freq_tables$sx_freq |>
  select(symptom_group, N_sx) |>
  group_by(symptom_group) |>
  summarise(N_sxgrp = sum(N_sx)) |>
  mutate(prop = round(100*calculate_prop(N_sxgrp, nrow(data)), 2)) |>
  filter(symptom_group == "None recorded") |>
  gt() |>
  gtsave(paste0(output_dir, "sxgrp_prop.rtf"))


# Tables - Cumulative % of Cancer Sites ---- 
cumulative_freq_site <-
  data |> 
  group_by(cancer_site_desc) |> 
  summarise(freq = n()) |> 
  arrange(desc(freq)) |> 
  mutate(prop = 100*calculate_prop(freq, sum(freq)),
         cumulative = format_numb(cumsum(prop)),
         prop = format_numb(prop)) |>
  left_join(cancer_codelist_data) |>
  tibble::rownames_to_column(var = "row")

cumulative_freq_site |> 
  gt() |> 
  gtsave(paste0(output_dir, "cumulative_freq_site.rtf"))

write.csv(cumulative_freq_site |> select(-row), paste0(output_dir, "cumulative_freq_site.csv"))


# Tables ----
## ICD-10, cancer and cancer groups ----
cancer_class <- 
  data |> 
  select(site_icd10_o2, cancer_site_desc, cancer_group) |>
  distinct() |>
  arrange(site_icd10_o2) |>
  mutate(icd10_3dig = substr(site_icd10_o2, 1, 3))

cancer_class_3dig <-
  cancer_class |>
  select(icd10=icd10_3dig, cancer_site_desc, cancer_group) |>
  distinct()

cancer_class_4dig <- cancer_class_3dig |> 
  group_by(icd10) |> 
  filter(n() > 1) |> 
  distinct(icd10) |>
  left_join(cancer_class, by=c("icd10"="icd10_3dig")) |>
  select(-icd10, icd10=site_icd10_o2)

final_cancer_class <- 
  cancer_class_3dig |> group_by(icd10) |> filter(n()==1) |> ungroup() |>
  bind_rows(cancer_class_4dig) |>
  arrange(icd10) |>
  distinct()

saveRDS(final_cancer_class, "./data/final_cancer_class.rds")

rm(cancer_class, cancer_class_3dig, cancer_class_4dig, final_cancer_class); gc()

## Symptom and Symptom groups ----
final_symptom_class <- 
  symptom_groups |>
  mutate(symptom_group = factor(symptom_group, c("Upper abdominal", "Lower abdominal",
                                                 "Breast Symptoms", "Central nervous system", "Lump/mass/lymph node", "Musculoskeletal",
                                                 "Respiratory", "Skin Lesion", "Ulceration", "Urological",
                                                 "Female specific", "Male specific", "Non-specific", "None recorded"))) |>
  arrange(symptom_group) |>
  select(symptom_group, symptom) |>
  rename(Symptom = symptom,
         "Symptom Group" = symptom_group)

saveRDS(final_symptom_class, "./data/final_symptom_class.RDS")


## Sample Characteristics ----
# information needed to produce table 1
table_info <-
  data |>
  inner_join(cohort |> select(patient_pseudoid, fiveyearageband)) |>
  select(cancer_site_desc, cancer_group, 
         sex, fiveyearageband) |>
  mutate(var = case_when(
    fiveyearageband %in% c("25-29", "30-34", "35-39",
                           "40-44", "45-49", "50-54", "55-59") ~ "<60",
    fiveyearageband %in% c("60-64", "65-69", "70-74") ~ "60-74",
    fiveyearageband %in% c("75-79", "80-84", "85-89", "90+") ~ "75+")) 

table1 <-
  table_info|>
  group_by(sex, var) |>
  summarise(n = n()) |>
  mutate(N = sum(n),
         cancer_site_desc = "Total",
         cancer_group = "Total",
         summary = str_glue({"{n} ({round(100*n/N, 1)}%)"})) |>
  ungroup() |>
  bind_rows(table_info |>
              group_by(cancer_site_desc, cancer_group, sex, var) |>
              summarise(n = n()) |> # numb per group
              mutate(N = sum(n),
                     summary = str_glue("{n} ({format_numb(100*n/N)}%)")) |>
              ungroup()) |>
  pivot_wider(id_cols = c(cancer_site_desc, cancer_group),
              values_from = c(summary, N),
              names_from = c(sex, var)) |>
  select(-c("N_1_<60", "N_1_60-74", "N_2_<60", "N_2_60-74")) |> # repeated columns
  rename(`N_1` = `N_1_75+`,
         `N_2` = `N_2_75+`) 

saveRDS(table1, file = "./data/table1_info.rds")

rm(table_info, table1)

## 'Numb' of symptoms per cancer site at different thresholds ----
## # information needed to produce table 2 and supplementary table 5
numb_sx_perc <-
  freq_tables$sx_site |>
  mutate(top1 = prop > 1 & symptom_group != "None recorded",
         top5 = prop > 5 & symptom_group != "None recorded",
         top10 = prop > 10 & symptom_group != "None recorded",
         top20 = prop > 20 & symptom_group != "None recorded",
         top50 = prop > 50 & symptom_group != "None recorded") |>
  group_by(cancer_site_desc, cancer_group, N_site) |>
  summarise(">1%" = sum(top1),
            ">5%" = sum(top5),
            ">10%" = sum(top10),
            ">20%" = sum(top20),
            ">50%" = sum(top50)) |>
  inner_join(freq_tables$sx_site |> 
               filter(symptom_group == "None recorded") |>
               group_by(cancer_site_desc, cancer_group, N_site) |>
               summarise(no_sx_prop = sum(prop)) |>
               mutate(prop = 100 - no_sx_prop) |>
               mutate(lower = 100*calculate_ci(prop=prop/100, N=N_site)$lower,
                      upper = 100*calculate_ci(prop=prop/100, N=N_site)$upper,
                      lower = ifelse(lower<0, 0, lower)) |>
               select(-c(no_sx_prop, N_site)) |>
               mutate(prop = format_numb(prop),
                      lower = format_numb(lower),
                      upper = format_numb(upper),
                      summary = str_glue("{prop} ({lower}, {upper})")) |>
               select(-c(prop, lower, upper)),
             by = c("cancer_site_desc"="cancer_site_desc", "cancer_group"="cancer_group")) |>
  arrange(cancer_site_desc) |>
  mutate(cancer_group = factor(cancer_group, levels = cancergrp_order)) |>
  rename(N = N_site) |>
  as_tibble() |>
  select(cancer_group, cancer_site_desc, N, `% with at least one recorded symptom`=summary, everything())

saveRDS(numb_sx_perc, file = "./data/numb_sx_perc.rds")

rm(numb_sx_perc)

## Proportions and CI for sx/site groups (htmap group) long format table ----
sx_site_grp_table <- freq_tables$sx_site_group |>
  mutate(lower = 100*calculate_ci(prop=prop/100, N=N_site)$lower,
         upper = 100*calculate_ci(prop=prop/100, N=N_site)$upper,
         lower = ifelse(lower<0, 0, lower)) |>
  mutate(prop = format_numb(prop),
         lower = format_numb(lower),
         upper = format_numb(upper), 
         summary = str_glue("{prop} ({lower}, {upper})"),
         symptom_group = factor(symptom_group, levels = sxgrp_order),
         cancer_group = factor(cancer_group, levels = cancergrp_order)) |>
  mutate(summary = ifelse(N_site=="N<20", "-", summary))  |>
  select(-c(n, prop, lower, upper))

saveRDS(sx_site_grp_table, "./data/sx_site_grp_table.rds")

rm(sx_site_grp_table)

## Proportions and CI for site/sx groups long format table (stratified by sex)----
site_sx_grp_table <- freq_tables$site_sx_group |>
  filter(N_sx>0) |>
  mutate(lower = 100*calculate_ci(prop=prop/100, N=N_sx)$lower,
         upper = 100*calculate_ci(prop=prop/100, N=N_sx)$upper,
         lower = ifelse(lower<0, 0, lower)) |>
  mutate(prop = format_numb(prop),
         lower = format_numb(lower),
         upper = format_numb(upper), 
         summary = str_glue("{prop} ({lower}, {upper})"),
         symptom_group = factor(symptom_group, levels = sxgrp_order),
         cancer_group = factor(cancer_group, levels = cancergrp_order),
         sex = ifelse(sex == 1, "Men", "Women")) |>
  mutate(summary = ifelse(N_sx=="N<20", "-", summary)) |>
  select(-c(n, prop, lower, upper))

saveRDS(site_sx_grp_table, "./data/site_sx_grp_table.rds")

rm(site_sx_grp_table)

## Proportions and CI for sx/site (heatmap) long format table ---- 
sx_site_table <- freq_tables$sx_site |>
  mutate(lower = 100*calculate_ci(prop=prop/100, N=N_site)$lower,
         upper = 100*calculate_ci(prop=prop/100, N=N_site)$upper,
         N_site = ifelse(N_site<20, "N<20", N_site)) |>
  mutate(
    lower = ifelse(lower<0, format_numb(0), format_numb(lower)), #lower bound to be 0 if neg numb
    upper = format_numb(upper),
    prop = format_numb(prop),
    summary = str_glue("{prop} ({lower}, {upper})"),
    `Cancer Site` = str_glue("{cancer_site_desc} (N = {N_site})")) |>
  mutate(symptom_group = factor(symptom_group, levels = sxgrp_order),
         summary = ifelse(N_site=="N<20", "-", summary))  |>
  select(-c(n, N_site, prop, lower, upper))

saveRDS(sx_site_table, "./data/sx_site_table.rds")
rm(sx_site_table)

## Proportions and CI for site/sx (barplot) long format table (stratified by sex) ----
site_sx_table <- 
  data_long |>
  group_by(sex, symptom, cancer_site_desc) |>
  summarise(n = sum(present)) |>
  mutate(N = sum(n), 
         prop = 100*n/N)|>
  mutate(lower = 100*calculate_ci(prop=prop/100, N=N)$lower,
         upper = 100*calculate_ci(prop=prop/100, N=N)$upper) |>
  mutate(lower = ifelse(lower<0, format_numb(0), format_numb(lower)),
         upper = format_numb(upper),
         prop = format_numb(prop),
         summary = str_glue("{prop} ({lower}, {upper})")) |>
  ungroup() |>
  inner_join(symptom_groups) |>
  inner_join(cancer_codelist_data |> select(-icd10_3dig, -icd10_4dig))

site_sx_table <-
  site_sx_table |>
  mutate(N = ifelse(N<20, "N<20", N)) |>
  arrange(symptom, cancer_site_desc) |>
  mutate(symptom_group = factor(symptom_group, levels = sxgrp_order),
         cancer_group = factor(cancer_group, levels = cancergrp_order)) |>
  mutate(symptom = factor(symptom, levels = unique(symptom[order(symptom_group)])),
         cancer_site_desc = factor(cancer_site_desc, levels = cancer_order),
         sex = ifelse(sex == 1, "Men", "Women"),
         summary = ifelse(N=="N<20", "-", summary)) |>
  select(-c(n, prop, lower, upper))

saveRDS(site_sx_table, "./data/site_sx_table.rds")

rm(site_sx_table)

## ICD-10, cancer site and cancer group classification ----
cancer_class <- data |> 
  group_by(cancer_site_desc) |> 
  summarise(freq = n()) |> 
  left_join(cancer_codelist_data) |>
  mutate(cancer_site_desc = factor(cancer_site_desc, levels = cancer_order),
         cancer_group = factor(cancer_group, levels = cancergrp_order)) |>
  arrange(cancer_group, cancer_site_desc)

cancer_class |>
  gt() |>
  gtsave(paste0(output_dir, "cancer_classification.rtf"))

write.csv(cancer_class, paste0(output_dir, "cancer_classification.csv"))

rm(cancer_class); gc()

# Data for Figures ----
sx_fig_name <- read.delim("../lookup_files/sx_figures_rename.txt")
## Heatmaps ----
{
  ### cancer site and sx level ----
  htmap_cancersite <- freq_tables$sx_site |>
    mutate(cancer_group = factor(cancer_group, levels = cancergrp_order)) |>
    distinct(cancer_group, cancer_site_desc) |>
    group_by(cancer_group) |> # 1 row per cancer group
    summarise(n = n(), # numb of cancer sites per group
              sites = paste(cancer_site_desc, collapse = ", ")) |>
    mutate(col = rep(c("#000000", "#7D7566"), length.out = n())) |>
    rowwise() |>
    mutate(colours = list(rep(col, n)))
  htmap_cancersite_col <- unlist(htmap_cancersite$colours)
  saveRDS(htmap_cancersite_col, file = "./data/htmap_cancersite_col.rds")
  
  htmap_sx <- freq_tables$sx_site |>
    mutate(symptom_group = factor(symptom_group, levels = sxgrp_order)) |>
    distinct(symptom_group, symptom) |>
    group_by(symptom_group) |> # 1 row per sx group
    summarise(n = n(), # numb of sx per group
              sites = paste(symptom, collapse = ", ")) |>
    mutate(col = rep(c("#000000", "#7D7566"), length.out = n())) |>
    rowwise() |>
    mutate(colours = list(rep(col, n)))
  htmap_sx_col <- rev(c(unlist(htmap_sx$colours), # rev so colours appear as need in htmap (top to bottom)
                        "#000000")) # black for 'avg no. symptom' label
  saveRDS(htmap_sx_col, file = "./data/htmap_sx_col.rds")
  
  
  # get average number of symptoms per cancer site 
  # exc. N/A and N/K (assumed to be no symptoms)
  htmap_df <-
    freq_tables$sx_site |>
    filter(cancer_group != 'N') |>
    mutate(prop = round(prop, 0)) |>
    bind_rows(data |>
                filter(`N/A` != 1 & `N/K` != 1) |> 
                group_by(cancer_site_desc) |>
                summarise(N_site = n(),
                          n = sum(numb_sx)) |>
                mutate(prop = round(calculate_prop(n, N_site), 1),
                       symptom = "Mean no. symptoms",
                       symptom_group = "summary")) |>
    mutate(symptom_recod = mapvalues(symptom,
                                     from = sx_fig_name$Original,
                                     to = sx_fig_name$New)) |> #renaming
    arrange(symptom_group) |>
    filter(N_site >= 20) |> # all cancer sites remain
    select(-c(n, N_site))
  
  saveRDS(htmap_df, file = "./data/htmap_df.rds")
  
  
  rm(list = ls(pattern = "^htmap"))
  
  
  ## cancer and symptom groups level ----
  htmap_grp_txt_col <- ifelse(round(freq_tables$sx_site_group$prop, 0) == 0, "#d3d3d3", "black")
  saveRDS(htmap_grp_txt_col, file = "./data/htmap_grp_txt_col.rds")
  
  
  htmap_grp_df <- freq_tables$sx_site_group |>
    filter(symptom_group != "summary") |>
    mutate(prop = round(prop, 0)) |>
    select(-c(n, N_site))
  saveRDS(htmap_grp_df, file = "./data/htmap_grp_df.rds")
  
  rm(list = ls(pattern = "^htmap")) 
}

## Barplots ----
sxgrp_fig_name <- read.delim("../lookup_files/sxgrp_figures_rename.txt")
sx_fig_name <- read.delim("../lookup_files/sx_figures_rename.txt")

### Site by symptom stratified by sex ----
barplt_df <-
  data_long |>
  group_by(sex, symptom, cancer_site_desc) |>
  summarise(n = sum(present)) |>
  mutate(N = sum(n), 
         prop = 100*calculate_prop(n, N))|>
  inner_join(symptom_groups) |>
  inner_join(cancer_codelist_data |> select(-icd10_3dig, -icd10_4dig)) |>
  ungroup() |>
  mutate(symptom_recod = mapvalues(symptom,
                                   from = sx_fig_name$Original,
                                   to = sx_fig_name$New),
         sex = ifelse(sex == 1, "Men", "Women"),
         symptom_group = mapvalues(symptom_group,
                                   from = sxgrp_fig_name$Old,
                                   to = sxgrp_fig_name$New)) |> #recoding
  arrange(cancer_site_desc)

barplt_df <-
  barplt_df |>
  filter(N >= 20)

saveRDS(barplt_df |> select(-c(n, N)), "./data/barplt_df.rds")

legend_grp <- # helper df for legend groups
  barplt_df |> 
  distinct(cancer_group, cancer_site_desc) |> 
  group_by(cancer_group) |> # 1 row per cancer group
  summarise(n = n(), # numb of cancer sites per group
            sites = paste(cancer_site_desc, collapse = ", ")) |> 
  mutate(to = cumsum(n), # indices that will be used in 'vector' to determine location of sites
         from = to - n + 1,
         palettes = c("red", "pink", "purple", "amber",
                      "indigo", "blue", "light-blue", "green",
                      "light-green", "teal", "red", "lime",
                      "brown", "brown", "orange", "grey"), # palette name per cancer group
         rev = c(rep(c(F,T), 6), F, F, F, T))
saveRDS(legend_grp, "./data/barplt_legend_grp.rds")

labels <- levels(barplt_df$cancer_site_desc)
labels <- setNames(labels, labels)
saveRDS(labels, "./data/barplt_labels.rds")

# create a palette for each cancer group using prespecified info from dataframe then merge colours into a vector
colors <- unlist(lapply(1:nrow(legend_grp), FUN = function(x) {
  c(pal_material(legend_grp$palettes[x], reverse = legend_grp$rev[x])(legend_grp$n[x]))})) 
colors <- setNames(colors, levels(barplt_df$cancer_site_desc))
saveRDS(colors, "./data/barplt_colors.rds")

### Cancer group by symptom group stratified by sex ----
barplt_grp_df <-  freq_tables$site_sx_group |> 
  mutate(sex = ifelse(sex == 1, "Men", "Women")) |>
  mutate(symptom_group = factor(symptom_group, levels = rev(sxgrp_order))) |>
  select(-c(n, N_sx))
saveRDS(barplt_grp_df, "./data/barplt_grp_df.rds")


# ########## END # ##########
