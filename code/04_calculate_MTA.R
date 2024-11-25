# # # # # # # # # # # # # # # #
# Purpose of this script is to calculate the MTA indicator for temporal
# period. This script assumes that the data has been preprocessed prior.
#
# Author:
# Martin Jung (IIASA) - 2024
# # # # # # # # # # # # # # # #

# Load Packages
library(sf)
library(dplyr)
library(tidyverse)
library(stringr)
library(assertthat)
library(purrr)
library(tidyr)
library(lubridate)

# Version - Update after each major correction.
version <- "01"

# Load custom functions
source("code/00_functions.R")

# path to the output
path_output <- "export/"
dir.create(path_output, showWarnings = FALSE)

# Create dataset with relevant files.
# Edit as needed
data <- data.frame(ifname = c("temporary_data/Speciesareas__Art17_habitats_distribution_2013_2018_EU.rds",
                            "temporary_data/Speciesareas__Art17_species_distribution_2013_2018_EU.rds",
                            "temporary_data/Speciesareas__EU_Art12_birds_distribution_2013_2018_with_sensitive_species.rds",
                   "temporary_data/Speciesareas__Art17_species_distribution_2013_2018_MS.rds",
                   "temporary_data/Speciesareas__Art17_habitats_distribution_2013_2018_MS.rds"),
                 scale = c("EU", "EU",  "EU", "MS", "MS"),
                 art = c("17","17","12", "17", "17"),
                 variable = c("habitats", "species", "species", "species", "habitats")
                 )
assertthat::assert_that(
  nrow(data)>0, all(file.exists(data$ifname)),
  "species" %in% data$variable,
  is.character(version)
)

#### Overall summary statistics ####
# Simply calculate some summary statistics here for
# the dashboard

# Number of unique species per group and dataset
df <- data |>
  dplyr::pull(ifname) |> map_dfr(readRDS) |>
  dplyr::filter(ctype == "combined") |>
  dplyr::group_by(name, category) |>
  dplyr::summarise(nr_features = dplyr::n_distinct(code),
                   mean_conservedarea = mean(conservedarea_km2))
df$category[is.na(df$category)] <- "Species"
write.csv(df,"export/Overall_statistics.csv",row.names = FALSE)

# Average conserved area and proportion per EU MS
df <- dplyr::bind_rows(
  data |> dplyr::filter(scale == "MS"),
  data |> dplyr::filter(art == 12)
) |> dplyr::pull(ifname) |> map_dfr(readRDS)
df$category[is.na(df$category)] <- "Species"

df <- df |>
  dplyr::filter(country != "All") |>
  dplyr::group_by(country, ctype, category) |>
  dplyr::summarise(nr_features = dplyr::n_distinct(code),
                   mean_conservedarea = mean(conservedarea_km2),
                   prop = (sum(conservedarea_km2) / sum(totalarea_km2)) )

df <- df |> dplyr::group_by(country, ctype, category) |> dplyr::summarise(
  nr_features = sum(nr_features),
  mean_conservedarea = mean(mean_conservedarea),
  prop = mean(prop)
)
write.csv(df,"export/Overall_statistics_MS.csv",row.names = FALSE)

#### Overall across EU - Calculate MTA for conservation areas ####
# We calculate the MTA across all areas and per time-period
# And for each dataset.

# Filter to EU only for now
biodiversity <- data |> dplyr::filter(scale == "EU")
scale <- "EU"

## -- First calculate overall with all datasets combined  -- ##
# Load all files of that scale
df <- biodiversity |>
  dplyr::pull(ifname) |> map_dfr(readRDS)

# Check for duplicates and if found add category first
if(anyDuplicated(df$code)>0){
  # Set missing category first
  df$category[is.na(df$category)] <- "Species"

  df$code <- paste0(stringr::str_sub(df$category,start = 0,end = 1),
                    df$code[which(duplicated(df$code))]
  )
}
# Set category to all
df$category <- "All"

# Reset for Birds directive the country as we are interested in overall patterns
df$country[df$name=="EU_Art12_birds_distribution_2013_2018_with_sensitive_species"] <- "All"

# Recategorize groups
df$year <- factor(df$year,
                  levels = c(2000, 2006, 2012, 2018, 2024),
                  labels = c("<2000", "2000-2006", "2006-2012","2012-2018", "2018-2024"))

# ---------------- #
## Calculate the amount and proportion in conservation areas
# First overall
overall <- df |>
  dplyr::filter(ctype == "combined", country == "All",
                year == "2018-2024") |>
  dplyr::group_by(category, ctype,code) |>
  dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2,na.rm = TRUE),
                 totalarea_km2 = sum(totalarea_km2,na.rm = TRUE)) |>
  dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
assertthat::assert_that(nrow(overall)>0,
                        all( between(overall$fullprop, 0,1) ) )

# Per year
peryear <- df |> dplyr::filter(country == "All") |>
  dplyr::group_by(category, ctype, year, code) |>
  dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                 totalarea_km2 = sum(totalarea_km2)) |>
  dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
assertthat::assert_that(all( between(peryear$fullprop, 0,1) ))

# Calculate targets based on the overall range
tr_flat <- calc_targets(current_range = overall$total_conservedarea_km2,
                        potential_range = overall$totalarea_km2,
                        option = "flat", default_target = 0.3) |>
  dplyr::mutate(code = overall$code)

# Calculate Jung et al. targets that minimize extinction risk
tr_exrisk <- calc_targets(data = overall,
                          current_range = overall$total_conservedarea_km2,
                          potential_range = overall$totalarea_km2,
                          option = "extinctrisk") |>
  dplyr::mutate(code = overall$code)

# Calculate log-linear targets
tr_loglinear <- calc_targets(data = overall,
                             current_range = overall$total_conservedarea_km2,
                             potential_range = overall$totalarea_km2,
                             option = "loglinear") |>
  dplyr::mutate(code = overall$code)

# Remove units
tr_flat <- tr_flat |> units::drop_units()
tr_exrisk <- tr_exrisk |> units::drop_units()
tr_loglinear <- tr_loglinear |> units::drop_units()

# Checks
assertthat::assert_that(
  !anyNA(tr_flat$target_absolute),
  all(tr_flat$code %in% overall$code),
  !anyNA(tr_loglinear$target_absolute),
  all(tr_loglinear$code %in% overall$code),
  !anyNA(tr_exrisk$target_absolute),
  all(tr_exrisk$code %in% overall$code)
)

# Make a copy but anonymizing the codes
write.csv(tr_loglinear |>
            # dplyr::slice(which(overall$ctype=="combined")) |>
            dplyr::select(option, target_relative),
          paste0("export/Overall_targetsloglinear_",version,".csv"),
                 row.names = FALSE)

# --- #
# Calculate the average MTA across time periods for the various targets
tr1 <- peryear |> units::drop_units() |>
  dplyr::left_join(tr_flat, by = "code") |>
  dplyr::group_by(category, ctype, option, year, code) |>
  # dplyr::reframe(mta = mta(target_relative,fullprop)) |>
  dplyr::reframe(mta = mta(target_absolute = target_absolute,
                           conserved_absolute = total_conservedarea_km2)) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,ctype, category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

tr2 <- peryear |> units::drop_units() |>
  dplyr::left_join(tr_loglinear) |> dplyr::ungroup() |>
  dplyr::group_by(category, ctype, option, year, code) |>
  # dplyr::reframe(mta = mta(target_relative = target_relative,proportion = fullprop)) |>
  dplyr::reframe(mta = mta(target_absolute = target_absolute,
                           conserved_absolute = total_conservedarea_km2)) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,ctype, category,year) |>
  dplyr::reframe(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

tr3 <- peryear |> units::drop_units() |>
  dplyr::left_join(tr_exrisk) |>
  dplyr::group_by(category,ctype, option, year, code) |>
  # dplyr::reframe(mta = mta(target_relative,fullprop)) |>
  dplyr::reframe(mta = mta(target_absolute = target_absolute,
                           conserved_absolute = total_conservedarea_km2)) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,ctype, category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup() |> as.data.frame()

# Combine
tr <- dplyr::bind_rows(
  #tr1,
  tr2
  #tr3
  ) |>
  # Add variables
  dplyr::mutate(variable = "All",
                scale = {{scale}},
                art = "All") |>
  dplyr::mutate(art = paste0("Art_",art)) |>
  tidyr::drop_na()

# Checks
assertthat::assert_that(nrow(tr)>0, all(dplyr::between(tr$y,0,1)))

# Write the output
ofname <- paste0(path_output, "MTA_loglinear_",scale,"_All_All_",version,".csv")
write.csv(tr, ofname, row.names = FALSE)

# ----------- #
## Furthermore calculate the overall MTA contribution per region and per country
# Here we use the targets derived from above, calculating their contribution of
# by geographic region.

# assertthat::assert_that(nrow(overall)>0,
#                         nrow(tr_loglinear)>0)
#
# # Summarize per country
# country <- df |>
#   dplyr::filter(country != "All", year == "2018-2024") |>
#   dplyr::group_by(category,ctype, country, code) |>
#   dplyr::reframe(country_conservedarea_km2 = sum(conservedarea_km2),
#                  totalarea_km2 = sum(totalarea_km2)) |>
#   dplyr::mutate(fullprop = (country_conservedarea_km2 / totalarea_km2) |> units::drop_units() ) |>
# units::drop_units()

# ----------- #
# Overall per dataset (Art 17, species, habitat, etc...)
for(i in 1:nrow(biodiversity) ){

  ifname = biodiversity$ifname[i]
  art = biodiversity$art[i]
  variable = biodiversity$variable[i]
  scale = biodiversity$scale[i]

  # For all
  # Load the folder
  df <- readRDS(ifname)

  # Add category
  if(!utils::hasName(df, "category")) df$category <- "species"

  # Recategorize groups
  df$year <- factor(df$year,
                    levels = c(2000, 2006, 2012, 2018, 2024),
                    labels = c("<2000", "2000-2006", "2006-2012","2012-2018", "2018-2024"))

  # Filter to all for Art 17
  if("All" %in% unique(df$country)){
    df <- df |> dplyr::filter(country == "All")
  }

  ## Calculate the amount and proportion in conservation areas
  # First overall
  overall <- df |>
    dplyr::group_by(category,ctype, code) |>
    dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                   totalarea_km2 = sum(totalarea_km2)) |>
    dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
  assertthat::assert_that(all( between(overall$fullprop, 0,1) ))

  # Per year
  peryear <- df |>
    dplyr::group_by(category,ctype, year, code) |>
    dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                   totalarea_km2 = sum(totalarea_km2)) |>
    dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
  assertthat::assert_that(all( between(overall$fullprop, 0,1) ))

  # Calculate targets based on the overall range
  tr_flat <- calc_targets(current_range = overall$total_conservedarea_km2,
                          potential_range = overall$totalarea_km2,
                          option = "flat", default_target = 0.3) |>
    dplyr::mutate(code = overall$code)

  # Calculate Jung et al. targets that minimize extinction risk
  tr_exrisk <- calc_targets(data = overall,
                            current_range = overall$total_conservedarea_km2,
                            potential_range = overall$totalarea_km2,
                            option = "extinctrisk") |>
    dplyr::mutate(code = overall$code)

  # Calculate rodrigues et al. targets
  tr_loglinear <- calc_targets(data = overall,
                               current_range = overall$total_conservedarea_km2,
                               potential_range = overall$totalarea_km2,
                               option = "loglinear") |>
    dplyr::mutate(code = overall$code)

  # Remove units
  tr_flat <- tr_flat |> units::drop_units()
  tr_exrisk <- tr_exrisk |> units::drop_units()
  tr_loglinear <- tr_loglinear |> units::drop_units()

  # --- #
  # Calculate the average MTA across time periods for the various targets
  tr1 <- peryear |> units::drop_units() |>
    dplyr::left_join(tr_flat, by = "code") |>
    dplyr::group_by(category, ctype, option, year, code) |>
    dplyr::reframe(mta = mta(target_absolute = target_absolute,
                             conserved_absolute = total_conservedarea_km2)) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,ctype, category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  tr2 <- peryear |> units::drop_units() |>
    dplyr::left_join(tr_loglinear) |>
    dplyr::group_by(category,ctype,option, year, code) |>
    # dplyr::reframe(mta = mta(target_relative,fullprop) ) |>
    dplyr::reframe(mta = mta(target_absolute = target_absolute,
                             conserved_absolute = total_conservedarea_km2)) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,ctype, category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  tr3 <- peryear |> units::drop_units() |>
    dplyr::left_join(tr_exrisk) |>
    dplyr::group_by(category,ctype, option, year, code) |>
    # dplyr::reframe(mta = mta(target_relative,fullprop) ) |>
    dplyr::reframe(mta = mta(target_absolute = target_absolute,
                             conserved_absolute = total_conservedarea_km2)) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,ctype,category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  # Combine
  tr <- dplyr::bind_rows(
    #tr1,
    tr2
    #tr3
    ) |>
    # Add variables
    dplyr::mutate(variable = {{variable}},
                  scale = {{scale}},
                  art = {{art}}) |>
    dplyr::mutate(art = paste0("Art_",art))

  # Write the output
  ofname <- paste0(path_output, "MTA_loglinear_",scale,"_Art",art,"_",variable,"_",version,".csv")
  write.csv(tr, ofname, row.names = FALSE)
}
stop("DONE with all EU wide MTA computations...")

#### EU MS - Calculate MTA based on EU MS reporting ####
# We calculate the MTA across all areas and per time-period
# And for each dataset.

# Filter to EU only for now
biodiversity <- dplyr::bind_rows(
  data |> dplyr::filter(scale == "MS"),
  data |> dplyr::filter(art == 12)
) |> dplyr::mutate(scale = "MS")
scale <- "MS"
assertthat::assert_that(nrow(biodiversity)>2)

## -- First calculate overall with all datasets combined  -- ##
# Load all files of that scale
df <- biodiversity |>
  dplyr::pull(ifname) |> map_dfr(readRDS)

# Check for duplicates and if found add category first
if(anyDuplicated(df$code)>0){
  # Set missing category first
  df$category[is.na(df$category)] <- "Species"

  df$code <- paste0(stringr::str_sub(df$category,start = 0,end = 1),
                    df$code[which(duplicated(df$code))]
  )
}
# Set category to all
df$category <- "All"

# Recategorize groups
df$year <- factor(df$year,
                  levels = c(2000, 2006, 2012, 2018, 2024),
                  labels = c("<2000", "2000-2006", "2006-2012","2012-2018", "2018-2024"))

## Calculate the amount and proportion in N2k sites
# First overall
overall <- df |> dplyr::group_by(category,ctype, code) |>
  dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                 totalarea_km2 = sum(totalarea_km2)) |>
  dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
assertthat::assert_that(all( between(overall$fullprop, 0,1) ))

# Per year
peryear <- df |> dplyr::group_by(category,ctype, year, code) |>
  dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                 totalarea_km2 = sum(totalarea_km2)) |>
  dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
assertthat::assert_that(all( between(overall$fullprop, 0,1) ))

# Calculate targets based on the overall range
tr_flat <- calc_targets(current_range = overall$total_conservedarea_km2,
                        potential_range = overall$totalarea_km2,
                        option = "flat", default_target = 0.3) |>
  dplyr::mutate(code = overall$code)

# Calculate Jung et al. targets that minimize extinction risk
tr_exrisk <- calc_targets(data = overall,
                          current_range = overall$total_conservedarea_km2,
                          potential_range = overall$totalarea_km2,
                          option = "extinctrisk") |>
  dplyr::mutate(code = overall$code)

# Calculate rodrigues et al. targets
tr_loglinear <- calc_targets(data = overall,
                             current_range = overall$total_conservedarea_km2,
                             potential_range = overall$totalarea_km2,
                             option = "loglinear") |>
  dplyr::mutate(code = overall$code)

# Remove units
tr_flat <- tr_flat |> units::drop_units()
tr_exrisk <- tr_exrisk |> units::drop_units()
tr_loglinear <- tr_loglinear |> units::drop_units()

# --- #
# Calculate the average MTA across time periods for the various targets
tr1 <- peryear |> units::drop_units() |>
  dplyr::left_join(tr_flat, by = "code") |>
  dplyr::group_by(category, ctype, option, year, code) |>
  # dplyr::reframe(
  #   mta = mta(target_relative,fullprop)
  # ) |>
  dplyr::reframe(mta = mta(target_absolute = target_absolute,
                           conserved_absolute = total_conservedarea_km2)) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,ctype, category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

tr2 <- peryear |> units::drop_units() |>
  dplyr::left_join(tr_loglinear) |>
  dplyr::group_by(category,ctype, option, year, code) |>
  # dplyr::reframe(mta = mta(target_relative,fullprop)) |>
  dplyr::reframe(mta = mta(target_absolute = target_absolute,
                           conserved_absolute = total_conservedarea_km2)) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,ctype, category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

tr3 <- peryear |> units::drop_units() |>
  dplyr::left_join(tr_exrisk) |>
  dplyr::group_by(category,ctype, option, year, code) |>
  # dplyr::reframe(mta = mta(target_relative,fullprop)) |>
  dplyr::reframe(mta = mta(target_absolute = target_absolute,
                           conserved_absolute = total_conservedarea_km2)) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,ctype,category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

# Combine
tr <- dplyr::bind_rows(
  #tr1,
  tr2
  #, tr3
  ) |>
  # Add variables
  dplyr::mutate(variable = "All",
                scale = {{scale}},
                art = "All") |>
  dplyr::mutate(art = paste0("Art_",art))

# Write the output
ofname <- paste0(path_output, "MTA_loglinear_",scale,"_All_All_",version,".csv")
write.csv(tr, ofname, row.names = FALSE)

# ----------- #
for(i in 1:nrow(biodiversity) ){

  ifname = biodiversity$ifname[i]
  art = biodiversity$art[i]
  variable = biodiversity$variable[i]
  scale = biodiversity$scale[i]

  # For all
  # Load the folder
  df <- readRDS(ifname)

  # Add category
  if(!utils::hasName(df, "category")) df$category <- "species"

  # Recategorize groups
  df$year <- factor(df$year,
                    levels = c(2000, 2006, 2012, 2018, 2024),
                    labels = c("<2000", "2000-2006", "2006-2012","2012-2018", "2018-2024"))

  ## Calculate the amount and proportion in conserved sites
  # First overall
  overall <- df |> dplyr::group_by(category,ctype,code) |>
    dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                   totalarea_km2 = sum(totalarea_km2)) |>
    dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
  assertthat::assert_that(all( between(overall$fullprop, 0,1) ))

  # Per year
  peryear <- df |> dplyr::group_by(category,ctype,year, code) |>
    dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                   totalarea_km2 = sum(totalarea_km2)) |>
    dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
  assertthat::assert_that(all( between(overall$fullprop, 0,1) ))

  # Calculate targets based on the overall range
  tr_flat <- calc_targets(current_range = overall$total_conservedarea_km2,
                          potential_range = overall$totalarea_km2,
                          option = "flat", default_target = 0.3) |>
    dplyr::mutate(code = overall$code)

  # Calculate Jung et al. targets that minimize extinction risk
  tr_exrisk <- calc_targets(data = overall,
                            current_range = overall$total_conservedarea_km2,
                            potential_range = overall$totalarea_km2,
                            option = "extinctrisk") |>
    dplyr::mutate(code = overall$code)

  # Calculate rodrigues et al. targets
  tr_loglinear <- calc_targets(data = overall,
                               current_range = overall$total_conservedarea_km2,
                               potential_range = overall$totalarea_km2,
                               option = "loglinear") |>
    dplyr::mutate(code = overall$code)

  # Remove units
  tr_flat <- tr_flat |> units::drop_units()
  tr_exrisk <- tr_exrisk |> units::drop_units()
  tr_loglinear <- tr_loglinear |> units::drop_units()

  # --- #
  # Calculate the average MTA across time periods for the various targets
  tr1 <- peryear |> units::drop_units() |>
    dplyr::left_join(tr_flat, by = "code") |>
    dplyr::group_by(category, ctype, option, year, code) |>
    # dplyr::reframe(
    #   mta = mta(target_relative,fullprop)
    # ) |>
    dplyr::reframe(mta = mta(target_absolute = target_absolute,
                             conserved_absolute = total_conservedarea_km2)) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,ctype,category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  tr2 <- peryear |> units::drop_units() |>
    dplyr::left_join(tr_loglinear) |>
    dplyr::group_by(category,ctype,option, year, code) |>
    # dplyr::reframe(mta = mta(target_relative,fullprop)) |>
    dplyr::reframe(mta = mta(target_absolute = target_absolute,
                             conserved_absolute = total_conservedarea_km2)) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,ctype,category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  tr3 <- peryear |> units::drop_units() |>
    dplyr::left_join(tr_exrisk) |>
    dplyr::group_by(category,ctype,option, year, code) |>
    # dplyr::reframe(mta = mta(target_relative,fullprop)) |>
    dplyr::reframe(mta = mta(target_absolute = target_absolute,
                             conserved_absolute = total_conservedarea_km2)) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,ctype,category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  # Combine
  tr <- dplyr::bind_rows(
    #tr1,
    tr2
    #, tr3
    ) |>
    # Add variables
    dplyr::mutate(variable = {{variable}},
                  scale = {{scale}},
                  art = {{art}}) |>
    dplyr::mutate(art = paste0("Art_",art))

  # Write the output
  ofname <- paste0(path_output, "MTA_loglinear_",scale,"_Art",art,"_",variable,"_",version,".csv")
  write.csv(tr, ofname, row.names = FALSE)
}
stop("DONE with all EU MS MTA computations...")
