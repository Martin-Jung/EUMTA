# # # # # # # # # # # # # # # #
# This scripts loads the Favourable Reference Range estimates from Article 12
# and Article 17 reporting and prepares them for the MTA target calculations.
#
# Author:
# Martin Jung (IIASA) - 2024-2026
# # # # # # # # # # # # # # # #

# Load Packages
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(tidyverse)
library(stringr)
library(assertthat)
library(foreach)
library(tidyr)
library(lubridate)
library(ibis.iSDM)

# Load custom functions
source("code/00_functions.R")

# Path to downloaded FRR data
# Source: https://sdi.eea.europa.eu/catalogue/srv/eng/catalog.search#/metadata/9f71b3e3-f8ec-442b-a2d5-c3c190605ac4
path_art17 <- "temporary_data/Article17_2020_dataset_csv/"
assertthat::assert_that(dir.exists(path_art17),
                        file.exists(file.path(path_art17, "Article17_2020_speciesEUassessment.csv")),
                        msg = "Article 17 FRR data folder does not exist. Please download the data first.")
# Art 12
path_art12 <- "temporary_data/Art 12/"
assert_that(dir.exists(path_art12),
            file.exists(file.path(path_art12, "Article12_2020_birdsEUpopulation.csv")),
            msg = "Article 12 FRR data folder does not exist. Please download the data first.")

# Get species list from preprocessed representation data
ll <- list.files("temporary_data", pattern = "Speciesareas__.*\\.rds", full.names = TRUE)

# --- #
# Get and consolidate species names and match to code
splist_art17 <- read.csv(paste0(path_art17, "/Article17_2020_species_check_list.csv"))
hblist_art17 <- read.csv(paste0(path_art17, "/Article17_2020_habitats_check_list.csv"))
# Get Art 12
splist_art12 <- read.csv(paste0(path_art12, "/Article12_2020_bird_check_list.csv"))

# Output path
path_output <- "temporary_data/Feature_targets_frr.rds"
if(file.exists(path_output)){
  cli::cli_alert_info("FRR targets already exist at {path_output}. Skipping processing.")
} else {
  cli::cli_alert_info("Processing FRR targets and saving to {path_output}...")
}

# -------------- #
#### Match species codes with FRR data ####

# Load Art 17 frr values for species and format
art_frr <- read.csv(paste0(path_art17, "/Article17_2020_speciesEUassessment.csv")) |>
  dplyr::select(country, region,
                assessment_speciescode, assessment_speciesname,
                range_surface_area,
                complementary_favourable_range, range_trend, complementary_favourable_range,
                complementary_favourable_range_q)

hab_frr <- read.csv(paste0(path_art17, "/Article17_2020_habitatsEUassessment.csv")) |>
  dplyr::select(country, region, habitatcode, habitatname,
                range_surface_area,
                complementary_favourable_range, range_trend, complementary_favourable_range,
                complementary_favourable_range_q)

