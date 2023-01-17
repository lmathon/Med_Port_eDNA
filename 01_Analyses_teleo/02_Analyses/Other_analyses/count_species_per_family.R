## ---------------------------------
## Count distinct species per family (Spygen teleo result)
##
## Author: Morgane Bruno
## Email: morgane.bruno@cefe.cnrs.fr
##
## Date Created: 2023-01-16
## ---------------------------------

#_libs
library(dplyr)
library(readr)
library(janitor)
library(taxize)
library(stringr)
library(vctrs)

#_main
## In port
taxo_teleo_port <- read.csv("01_Analyses_teleo/00_data/Teleo_SC21250_résultats_campagnes 1 et 2.csv", sep=";") %>%
  dplyr::as_tibble() %>%
  dplyr::select(project, X, X.1, X.2) %>%
  janitor::row_to_names(dat = ., row_number = 4, remove_rows_above = T, remove_row = T) %>%
  dplyr::mutate(scientific_name = gsub(pattern = ' ', replacement = '_', x = scientific_name))

### Number of species per family
port_count_family <- readr::read_csv("01_Analyses_teleo/00_data/matrice_teleo_port.csv", show_col_types = FALSE) %>%
  dplyr::select(...1) %>%
  dplyr::rename_with(.cols = 1, .fn = ~"scientific_name") %>%
  dplyr::mutate(scientific_name = gsub(pattern = ' ', replacement = '_', x = scientific_name)) %>%
  dplyr::left_join(x = ., y = taxo_teleo_port, by = 'scientific_name') %>%
  dplyr::filter(stringr::str_count(string = scientific_name, pattern = "_") != 0) %>% # Remove Taxon if is family_name
  dplyr::group_by(family) %>%
  dplyr::distinct() %>% 
  dplyr::summarise(species_count = n()) %>%
  dplyr::mutate(sampling = "Port")

## Outside Port, Biohut
taxo_teleo_biohut <- readr::read_csv("01_Analyses_teleo/00_data/biodiv_milieu_naturel.csv", show_col_types = FALSE) %>%
  dplyr::select(...1) %>%
  dplyr::rename_with(.cols = 1, .fn = ~"scientific_name") %>%
  dplyr::mutate(scientific_name = gsub(pattern = ' ', replacement = '_', x = scientific_name)) %>%
  dplyr::left_join(x = ., y = taxo_teleo_port, by = 'scientific_name')

### Request NCBI Taxonomy for Taxon whose family == NA
taxon_to_query <- taxo_teleo_biohut %>% 
  dplyr::filter(is.na(family)) %>% 
  dplyr::mutate(taxon = dplyr::case_when(
    # keep only Genus, if Taxon == Genus sp. OR Taxon == multiple species
    grepl(pattern = "_sp.$", x = scientific_name) ~ stringr::word(string = scientific_name, start = 1, end = 1, sep = "_"),
    stringr::str_count(string = scientific_name, pattern = "_") > 1 ~ stringr::word(string = scientific_name, start = 1, end = 1, sep = "_"),
    TRUE ~ sub(pattern = "_", replacement = " ", x = scientific_name)
  ))
taxon_to_query[taxon_to_query$taxon == "C.",]$taxon <- "Hirundichthys speculiger" # replace C. from C._heterurus_H._speculiger by Cheilopogon
res_query <- taxon_to_query %>%
  dplyr::mutate(chunk = vctrs::vec_rep_each(c(1,2,3,4), 50)[1:155]) %>%
  dplyr::group_split(chunk) %>%
  lapply(., function(subset){
    family <- taxize::tax_name(sci = subset$taxon, get = "family", db = "ncbi")$family
    return(family)
  })
taxon_to_query <- taxon_to_query %>% 
  dplyr::mutate(family = res_query %>% unlist()) %>%
  dplyr::select(-taxon)

#### Number of species per family
biohut_count_family <- taxo_teleo_biohut %>%
  dplyr::filter(!is.na(family)) %>%
  rbind(., taxon_to_query) %>%
  dplyr::filter(stringr::str_count(string = scientific_name, pattern = "_") != 0) %>% # Remove Taxon if is family_name
  dplyr::group_by(family) %>%
  dplyr::distinct() %>% 
  dplyr::summarise(species_count = n()) %>%
  dplyr::mutate(sampling = "Biohut")


write.csv(x = rbind(port_count_family, biohut_count_family), file = "20230116_port_biohut_count_species_per_family.csv")
