## ---------------------------
##
## Summary table - species list
##
## Author: Morgane Bruno
## Email: morgane.bruno@cefe.cnrs.fr
##
## Date Created: 2023-01-20
## ----------------------------

#_libs
library(readxl)
library(dplyr)
library(tibble)
library(gt)

#_input
expert_file <- '01_Analyses_teleo/00_data/ListePoissonsExpert.csv'
cate_file <- '01_Analyses_teleo/03_Outputs/Species_list_per_category.csv'
port_file <- '01_Analyses_teleo/00_data/matrice_teleo_port.csv'
reserve_file <- '01_Analyses_teleo/03_Outputs/Species_list_reserve.csv'
lockdown_file <- '01_Analyses_teleo/03_Outputs/Species_list_lockdown.csv'
trait_file <- '01_Analyses_teleo/00_data/Functional_data_corrected_20220124.csv'

#_main

## 01-Load Data
sp_expert_port <- read.csv(file = expert_file, header = FALSE, col.names = 'Species') %>%
  dplyr::mutate(Species = gsub(pattern = '_', replacement = ' ', x = Species))
sp_only_cate <- read.csv(file = cate_file)
sp_port <- read.csv(file = port_file, row.names = 1) %>%
  tibble::rownames_to_column(var = 'Species') %>% 
  dplyr::select(Species) %>%
  dplyr::mutate(Species = gsub(pattern = '_', replacement = ' ', x = Species))
sp_reserve <- read.csv(file = reserve_file, header = TRUE) %>%
  dplyr::mutate(Species = gsub(pattern = '_', replacement = ' ', x = Species))
sp_lockdown <- read.csv(file = lockdown_file, header = TRUE) %>%
  dplyr::mutate(Species = gsub(pattern = '_', replacement = ' ', x = Species))
traits <- read.csv(file = trait_file, header = TRUE) %>%
  dplyr::select(Species, all_commercial_level, IUCN_Red_List_Category, CryptoBenthic)

## 02-Create summary table
sp_notport <- sp_reserve %>% dplyr::filter(!Species %in% sp_port$Species)
summ <- sp_port %>%
  # According to expert: species present in ports or not
  dplyr::mutate(Expert = as.integer(Species %in% sp_expert_port$Species)) %>%
  # Add species not found in port (reserve)
  dplyr::bind_rows(sp_notport) %>%
  # Add trait columns: commercial, iucn, cryptic
  dplyr::left_join(x = ., y = traits, by = 'Species') %>%
  dplyr::mutate(
    # Convert iucn category to threatened or not
    Threatened = dplyr::case_when(
      IUCN_Red_List_Category %in% c('VU', 'EN', 'CR') ~ 1, # Vulnerable, Endangered, Critically Endangered
      IUCN_Red_List_Category %in% c('DD', 'LC', 'NT') ~ 0  # Data Deficient, Least Concern, Near Threatened
    ),
    Port = dplyr::case_when(
      Species %in% sp_only_cate$Port_only ~ 2, # only in reserve
      Species %in% sp_port$Species ~ 1, #present in port
      TRUE ~ 0 #missing in reserve
    ),
    Reserve = dplyr::case_when(
      Species %in% sp_only_cate$Reserve_only ~ 2, # only in reserve
      Species %in% sp_reserve$Species ~ 1, #present in port
      TRUE ~ 0 #missing in reserve
    ),
    Lockdown = dplyr::case_when(
      Species %in% sp_only_cate$Lockdown_only ~ 2, # only in reserve
      Species %in% sp_lockdown$Species ~ 1, #present in port
      Species %in% sp_only_cate$Port_only ~ NA_integer_, # no observations in ports during lockdown
      TRUE ~ 0 #missing in reserve
    )
  ) %>%
  dplyr::select(-IUCN_Red_List_Category) %>%
  dplyr::arrange(desc(Expert))

## 03-Format output
colnames(summ) <- c('Species', 'In port (according to expert)', 'Commercial', 'Cryptic', 'Threatened', 'In port', 'In reserve', 'During lockdown')
summ <- summ %>% dplyr::select(c('Species', 'Cryptic', 'Commercial', 'Threatened', 'In port (according to expert)', 'In port', 'In reserve', 'During lockdown'))
write.csv(x = summ, file = '01_Analyses_teleo/03_Outputs/summary_table_traits_species.csv', row.names = FALSE)