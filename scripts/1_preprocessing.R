## --------------------------------
## Script name: 1_preprocessing.R
##
## Author: Nadine Zakkak
##
## Purpose: Preprocessing NCDA data
##
## ---------------------------------
# Loading libraries ----------------------
pacman::p_load(readr, plyr, dplyr, ggplot2, tidyr, fastDummies, gt, stringr, forcats)

# Set up working space and data -----------------
options(width = 80)
# ncda <- read.csv("path/to/data", stringsAsFactors = F)
colnames(ncda) <- tolower(colnames(ncda))
ncda <- ncda |> rename_with(.fn = ~gsub("ncda_", "", .x), cols=starts_with("ncda_"))
ncda_orig <- ncda
ncda <- ncda |> 
  select(tumour_pseudoid, patient_pseudoid, pseudo_gp_code,
         sex, fiveyearageband, quintile_2019,ca_reg_ethnicity_group,
         site_icd10_o2, morph_coded, morph_icd10_o2, behaviour_icd10_o2, behaviour_coded,
         grade, stage_best, stage_best_system,
         presentsymptom, presentsymptom_desc, presentsignortest, presentsignortest_desc,
         investigations, investigations_desc,
         typereferral, typereferral_desc,
         diagnostic_interval, patient_interval, pci, pci_alt, ref_to_spec_cons,
         delaypresentation, delayafterreferral,
         chrl_tot_27_03, chrl_tot_78_06,
         screening_participation, screening_participationd,
         screen_detected, screen_detected_desc,
         comorbidity, comorbidity_desc,
         placepresentation, placepresentation_desc,
         consultations, multiconsults, multiconsultsd, consultations_unknown) |> 
  dplyr::rename(stage = stage_best,
                symptom = presentsymptom,
                symptom_desc = presentsymptom_desc)

patients_numb <- data.frame(description = "original", 
                            patients = length(unique(ncda$patient_pseudoid)), 
                            tumours = length(unique(ncda$tumour_pseudoid)),
                            stringsAsFactors = FALSE) #to save cohort numbers

output_dir <- "./results/DI/"

sx_recode <- read.delim("./lookup_files/sx_recode.txt", sep = " ", header = TRUE)
sx_recode_desc <- read.delim("./lookup_files/sx_recode_desc.txt", header=TRUE)
symptom_groups <- read.delim("./lookup_files/symptom_groups.txt", sep=" ")
sx_fig_name <- read.delim("./lookup_files/sx_figures_rename.txt")

cancer_groups <- read.delim("./lookup_files/cancer_site_groups.txt") # THIS CODELIST WAS BASED ON THE DATA AVAILABLE - DOESN'T CONTAIN ALL ICD-10 CODES **

# Convert to 10 year age bands ---- 
ncda <-
  ncda |> 
  mutate(
    tenyearageband = case_when(
      fiveyearageband %in% c("0-4", "5-9", "10-14", "15-19", "20-24") ~ "<25",
      fiveyearageband %in% c("25-29") ~ "25-29",
      fiveyearageband %in% c("30-34", "35-39") ~ "30-39",
      fiveyearageband %in% c("40-44", "45-49") ~ "40-49",
      fiveyearageband %in% c("50-54", "55-59") ~ "50-59",
      fiveyearageband %in% c("60-64", "65-69") ~ "60-69",
      fiveyearageband %in% c("70-74", "75-79") ~ "70-79",
      fiveyearageband %in% c("80-84", "85-89") ~ "80-89",
      fiveyearageband == "90+" ~ "90+"
    ))


# Recode stage ----
ncda$stage_recode <- ncda$stage
ncda$stage_recode[which(ncda$stage_recode %in% c("?", "U" , "", "A", "B", "C"))] <- -99 
ncda$stage_recode <- parse_number(ncda$stage_recode) #only extract number


# Incl./Excl. steps ----
ncda_backup <- ncda
## Excl screen-detected cancers
ncda <- ncda |> filter(screen_detected != 1 & typereferral != 5)

patients_numb <- rbind(patients_numb,
                       c("Not via screening", 
                         length(unique(ncda$patient_pseudoid)),
                         length(unique(ncda$tumour_pseudoid))))

## Only patients aged 25+
ncda <- ncda |> filter(tenyearageband != "<25") 

patients_numb <- rbind(patients_numb,
                       c("Age >=25", 
                         length(unique(ncda$patient_pseudoid)),
                         length(unique(ncda$tumour_pseudoid))))

## Diagnostic Interval (0,730) days or missing
ncda <- ncda |> filter(is.na(diagnostic_interval) | between(diagnostic_interval, 0, 730)) 

patients_numb <- rbind(patients_numb,
                       c("diagnostic interval (0, 730) days or missing",
                         length(unique(ncda$patient_pseudoid)),
                         length(unique(ncda$tumour_pseudoid))))

# Cancer sites Description 'super-groups' ---------- 
ncda_backup <- ncda
#C649 is not in the codelist, change to C64
setdiff(unique(ncda$site_icd10_o2),  cancer_groups$icd10_4dig)
ncda$site_icd10_o2[ncda$site_icd10_o2=="C649"] <- "C64"
ncda <-
  ncda |>
  inner_join(cancer_groups, by = c("site_icd10_o2"="icd10_4dig"))

# Incl./Excl. cont. steps ----
## Excl. Discordant sex - site
ncda_backup <- ncda
ncda <-
  ncda |>
  filter(!((sex == 1 & cancer_group == "Gynaecological")
           | (sex == 2 & cancer_group == "Prostate and other male organs")))

patients_numb <- rbind(patients_numb,
                       c("concordant sex-site", 
                         length(unique(ncda$patient_pseudoid)),
                         length(unique(ncda$tumour_pseudoid))))

# Recode symptom ----
ncda_backup <- ncda
# each sx on separate row
ncda_sx <- separate_rows(ncda, symptom, sep=" ")
#recoding
ncda_sx$symptom_recod <- as.numeric(mapvalues(ncda_sx$symptom,
                                              from = sx_recode$sx_old,
                                              to = sx_recode$sx_new))
# add description of symptom (per row)
ncda_sx <- ncda_sx |> 
  inner_join(sx_recode_desc |> dplyr::rename(symptom_recod_desc = sx_desc),
             by=c("symptom_recod"="sx")) 

ncda_sx <- ncda_sx |> 
  select(patient_pseudoid, tumour_pseudoid, pseudo_gp_code, sex,
         symptom_recod, symptom_recod_desc,
         stage_recode, site_icd10_o2, cancer_site_desc) 

ncda_sx <- unique(ncda_sx)

ncda_sx <-
  ncda_sx |>
  inner_join(symptom_groups, by = c("symptom_recod_desc"="symptom"))

# Incl./Excl. cont. steps ----
## Excl. Discordant sex - symptom 
discord_patid <- 
  ncda_sx |> 
  filter((sex == 1 & symptom_group == "Female specific")
         | (sex == 2 & symptom_group == "Male specific")) |>
  select(patient_pseudoid) |> pull()

ncda_sx <-
  ncda_sx |>
  filter(!(patient_pseudoid %in% discord_patid)) 

ncda <-
  ncda |>
  filter(!(patient_pseudoid %in% discord_patid))

patients_numb <- rbind(patients_numb,
                       c("concordant sex-symptom", 
                         length(unique(ncda$patient_pseudoid)),
                         length(unique(ncda$tumour_pseudoid))))

rm(sx_recode, sx_recode_desc)

## wide format
ncda_sx_wide <- dummy_cols(ncda_sx, 
                           select_columns=c("symptom_recod_desc"), 
                           remove_selected_columns = T)
ncda_sx_wide <- 
  ncda_sx_wide |> 
  group_by(patient_pseudoid, tumour_pseudoid, sex) |>
  dplyr::summarise(across(starts_with("symptom_recod_desc"), ~max(.x))) |>
  dplyr::mutate(numb_sx = rowSums(across(starts_with("symptom_recod_desc")))) |> 
  rename_with(.fn = ~gsub("symptom_recod_desc_", "symptom_", .x), 
              cols = starts_with("symptom_recod_desc_")) |>
  ungroup()

## long format
symptoms_long <- 
  ncda_sx_wide |> 
  pivot_longer(starts_with("symptom_"), 
               names_prefix = "symptom_",
               names_to = "symptom", 
               values_to = "present")

# multiple tumours recorded ----
# Choose random: a- highest stage
#                b- no highest stage --> randomly select
ncda_sx_wide <-
  ncda_sx_wide |> 
  left_join(ncda |> select(patient_pseudoid, tumour_pseudoid, 
                           site_icd10_o2 ,cancer_site_desc, stage_recode)) |> 
  as_tibble()

set.seed(154416)
single_tumour <-
  ncda_sx_wide |>
  group_by(patient_pseudoid) |>
  slice(which.max(rank(stage_recode, ties.method = "random", na.last = F))) |>
  as_tibble()

patients_numb <- rbind(patients_numb,
                       c("single tumour/patient", 
                         length(unique(single_tumour$patient_pseudoid)),
                         length(unique(single_tumour$tumour_pseudoid))))

# Add cancer site super groups ----
single_tumour <-
  single_tumour |>
  rename(orig_cancer_site_desc = cancer_site_desc) |>
  inner_join(cancer_groups, by = c("site_icd10_o2"="icd10_4dig"))

# Convert to long dataset ----
single_tumour_long <- 
  single_tumour |> 
  pivot_longer(cols = starts_with("symptom_"), names_prefix = "symptom_",
               names_to = "symptom", values_to = "present")

# Add symptom super groups ----
single_tumour_long <-
  single_tumour_long |>
  inner_join(symptom_groups, by = c("symptom"="symptom"))

# Rename columns in wide data set -----
single_tumour <- single_tumour |>
  rename_with(.fn = ~gsub("symptom_", "", .x), 
              cols = starts_with("symptom_"))

# Patients characteristics ----
cohort <-
  single_tumour |> 
  select(tumour_pseudoid, patient_pseudoid) |> 
  left_join(ncda, by = c('tumour_pseudoid' = 'tumour_pseudoid', "patient_pseudoid" = "patient_pseudoid"))

# Patient flowchart file ----
patients_numb <- patients_numb |>
  mutate(patients_diff = as.numeric(as.character(patients)) - lag(as.numeric(as.character(patients))),
         tumours_diff = as.numeric(as.character(tumours)) - lag(as.numeric(as.character(tumours)))) |>
  select(description, patients, patients_diff, tumours, tumours_diff)


# saveRDS(single_tumour, "/path/to/data/single_tumour_wide.RDS")
# saveRDS(single_tumour_long, "/path/to/data/single_tumour_long.RDS")
# saveRDS(cohort, "/path/to/data/cohort.RDS")
# write.csv(patients_numb, "/path/to/data/cohort_flowchart.csv", row.names = F)



