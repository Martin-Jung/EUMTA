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

# Load custom functions
source("code/00_functions.R")

# path to the output
path_output <- "export/"
dir.create(path_output, showWarnings = FALSE)

# Checks
assertthat::assert_that(
  file.exists("temporary_data/N2k.gpkg"),
  file.exists("temporary_data/Art17_2018.gpkg")
)

# Create dataset with relevant files.
# Edit as needed
data <- data.frame(ifname = c("temporary_data/Speciesareas__Art17_habitats_distribution_2013_2018_EU.rds",
                            "temporary_data/Speciesareas__Art17_habitats_distribution_2013_2018_MS.rds",
                            "temporary_data/Speciesareas__Art17_species_distribution_2013_2018_EU.rds",
                            "temporary_data/Speciesareas__Art17_species_distribution_2013_2018_MS.rds",
                            "temporary_data/Speciesareas__EU_Art12_birds_distribution_2013_2018_with_sensitive_species.rds"),
                 scale = c("EU", "MS", "EU", "MS", "EU"),
                 art = c("17","17","17","17","12"), variable = c("habitats", "habitats", "species", "species","species")
                 )
assertthat::assert_that(
  nrow(data)>0, all(file.exists(data$ifname)),
  "species" %in% data$variable
)

#### Overall across EU - Calculate MTA based on N2k areas ####
# We calculate the MTA across all areas and per time-period
# And for each dataset.

# Filter to EU only for now
biodiversity <- data |> dplyr::filter(scale == "EU")
scale <- "EU"

## -- First calculate overall with all datasets combined  -- ##
# Load all files of that scale
df <- biodiversity |>
  dplyr::pull(ifname) |> map_dfr(readRDS) |>
  # Recalculate to km2
  dplyr::mutate(conservedarea_km2 = units::set_units(conservedarea_m2, km^2),
                totalarea_km2 = units::set_units(totalarea_m2, km^2))

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
overall <- df |> dplyr::group_by(category,code) |>
  dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                 totalarea_km2 = sum(totalarea_km2)) |>
  dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
assertthat::assert_that(all( between(overall$fullprop, 0,1) ),
                        !any(duplicated(overall$code)))

# Per year
peryear <- df |> dplyr::group_by(category,year, code) |>
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
  dplyr::group_by(category, option, year, code) |>
  dplyr::reframe(
    mta = mean(
      1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
    )
  ) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

tr2 <- peryear |> units::drop_units() |>
  dplyr::left_join(tr_loglinear) |>
  dplyr::group_by(category,option, year, code) |>
  dplyr::reframe(mta = mean(
    1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
  )) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

tr3 <- peryear |> units::drop_units() |>
  dplyr::left_join(tr_exrisk) |>
  dplyr::group_by(category,option, year, code) |>
  dplyr::reframe(mta = mean(
    1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
  )) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

# Combine
tr <- dplyr::bind_rows(tr1, tr2, tr3) |>
  # Add variables
  dplyr::mutate(variable = "All",
                scale = {{scale}},
                art = "All") |>
  dplyr::mutate(art = paste0("Art_",art))

# Write the output
ofname <- paste0(path_output, "MTA_alltargets_",scale,"_All_All.csv")
write.csv(tr, ofname, row.names = FALSE)

# ----------- #
for(i in 1:nrow(biodiversity) ){

  ifname = biodiversity$ifname[i]
  art = biodiversity$art[i]
  variable = biodiversity$variable[i]
  scale = biodiversity$scale[i]

  # Skip
  if(calculated_all && art == "all") next()

  # For all
  # Load the folder
  df <- readRDS(ifname) |>
    # Recalculate to km2
    dplyr::mutate(conservedarea_km2 = units::set_units(conservedarea_m2, km^2),
                  totalarea_km2 = units::set_units(totalarea_m2, km^2))

  # Add category
  if(!utils::hasName(df, "category")) df$category <- "species"

  # Recategorize groups
  df$year <- factor(df$year,
                    levels = c(2000, 2006, 2012, 2018, 2024),
                    labels = c("<2000", "2000-2006", "2006-2012","2012-2018", "2018-2024"))

  ## Calculate the amount and proportion in N2k sites
  # First overall
  overall <- df |> dplyr::group_by(category,code) |>
    dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                   totalarea_km2 = sum(totalarea_km2)) |>
    dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
  assertthat::assert_that(all( between(overall$fullprop, 0,1) ),
                          !any(duplicated(overall$code)))

  # Per year
  peryear <- df |> dplyr::group_by(category,year, code) |>
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
    dplyr::group_by(category, option, year, code) |>
    dplyr::reframe(
      mta = mean(
        1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
        )
    ) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  tr2 <- peryear |> units::drop_units() |>
    dplyr::left_join(tr_loglinear) |>
    dplyr::group_by(category,option, year, code) |>
    dplyr::reframe(mta = mean(
                     1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
                   )) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  tr3 <- peryear |> units::drop_units() |>
    dplyr::left_join(tr_exrisk) |>
    dplyr::group_by(category,option, year, code) |>
    dplyr::reframe(mta = mean(
                     1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
                   )) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  # Combine
  tr <- dplyr::bind_rows(tr1, tr2, tr3) |>
    # Add variables
    dplyr::mutate(variable = {{variable}},
                  scale = {{scale}},
                  art = {{art}}) |>
    dplyr::mutate(art = paste0("Art_",art))

  # Write the output
  ofname <- paste0(path_output, "MTA_alltargets_",scale,"_Art",art,"_",variable,".csv")
  write.csv(tr, ofname, row.names = FALSE)
}
stop("DONE with all MTA computations...")


#### EU MS - Calculate MTA based on N2k areas ####
# We calculate the MTA across all areas and per time-period
# And for each dataset.

# Filter to EU only for now
biodiversity <- dplyr::bind_rows(
  data |> dplyr::filter(scale == "MS"),
  data |> dplyr::filter(art == 12)
) |> dplyr::mutate(scale = "MS")
scale <- "MS"

## -- First calculate overall with all datasets combined  -- ##
# Load all files of that scale
df <- biodiversity |>
  dplyr::pull(ifname) |> map_dfr(readRDS) |>
  # Recalculate to km2
  dplyr::mutate(conservedarea_km2 = units::set_units(conservedarea_m2, km^2),
                totalarea_km2 = units::set_units(totalarea_m2, km^2))

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
overall <- df |> dplyr::group_by(category,code) |>
  dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                 totalarea_km2 = sum(totalarea_km2)) |>
  dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
assertthat::assert_that(all( between(overall$fullprop, 0,1) ),
                        !any(duplicated(overall$code)))

# Per year
peryear <- df |> dplyr::group_by(category,year, code) |>
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
  dplyr::group_by(category, option, year, code) |>
  dplyr::reframe(
    mta = mean(
      1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
    )
  ) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

tr2 <- peryear |> units::drop_units() |>
  dplyr::left_join(tr_loglinear) |>
  dplyr::group_by(category,option, year, code) |>
  dplyr::reframe(mta = mean(
    1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
  )) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

tr3 <- peryear |> units::drop_units() |>
  dplyr::left_join(tr_exrisk) |>
  dplyr::group_by(category,option, year, code) |>
  dplyr::reframe(mta = mean(
    1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
  )) |>
  # Ungroup and then average per year
  dplyr::ungroup() |>
  dplyr::group_by(option,category,year) |>
  dplyr::summarise(y = mean( mta ),
                   sd = sd(mta)
  ) |> dplyr::ungroup()

# Combine
tr <- dplyr::bind_rows(tr1, tr2, tr3) |>
  # Add variables
  dplyr::mutate(variable = "All",
                scale = {{scale}},
                art = "All") |>
  dplyr::mutate(art = paste0("Art_",art))

# Write the output
ofname <- paste0(path_output, "MTA_alltargets_",scale,"_All_All.csv")
write.csv(tr, ofname, row.names = FALSE)

# ----------- #
for(i in 1:nrow(biodiversity) ){

  ifname = biodiversity$ifname[i]
  art = biodiversity$art[i]
  variable = biodiversity$variable[i]
  scale = biodiversity$scale[i]

  # Skip
  if(calculated_all && art == "all") next()

  # For all
  # Load the folder
  df <- readRDS(ifname) |>
    # Recalculate to km2
    dplyr::mutate(conservedarea_km2 = units::set_units(conservedarea_m2, km^2),
                  totalarea_km2 = units::set_units(totalarea_m2, km^2))

  # Add category
  if(!utils::hasName(df, "category")) df$category <- "species"

  # Recategorize groups
  df$year <- factor(df$year,
                    levels = c(2000, 2006, 2012, 2018, 2024),
                    labels = c("<2000", "2000-2006", "2006-2012","2012-2018", "2018-2024"))

  ## Calculate the amount and proportion in N2k sites
  # First overall
  overall <- df |> dplyr::group_by(category,code) |>
    dplyr::reframe(total_conservedarea_km2 = sum(conservedarea_km2),
                   totalarea_km2 = sum(totalarea_km2)) |>
    dplyr::mutate(fullprop = (total_conservedarea_km2 / totalarea_km2) |> units::drop_units() )
  assertthat::assert_that(all( between(overall$fullprop, 0,1) ),
                          !any(duplicated(overall$code)))

  # Per year
  peryear <- df |> dplyr::group_by(category,year, code) |>
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
    dplyr::group_by(category, option, year, code) |>
    dplyr::reframe(
      mta = mean(
        1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
      )
    ) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  tr2 <- peryear |> units::drop_units() |>
    dplyr::left_join(tr_loglinear) |>
    dplyr::group_by(category,option, year, code) |>
    dplyr::reframe(mta = mean(
      1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
    )) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  tr3 <- peryear |> units::drop_units() |>
    dplyr::left_join(tr_exrisk) |>
    dplyr::group_by(category,option, year, code) |>
    dplyr::reframe(mta = mean(
      1 - (pmin(1, pmax(0, ( target_relative - fullprop ) )) / target_relative )
    )) |>
    # Ungroup and then average per year
    dplyr::ungroup() |>
    dplyr::group_by(option,category,year) |>
    dplyr::summarise(y = mean( mta ),
                     sd = sd(mta)
    ) |> dplyr::ungroup()

  # Combine
  tr <- dplyr::bind_rows(tr1, tr2, tr3) |>
    # Add variables
    dplyr::mutate(variable = {{variable}},
                  scale = {{scale}},
                  art = {{art}}) |>
    dplyr::mutate(art = paste0("Art_",art))

  # Write the output
  ofname <- paste0(path_output, "MTA_alltargets_",scale,"_Art",art,"_",variable,".csv")
  write.csv(tr, ofname, row.names = FALSE)
}
stop("DONE with all MTA computations...")
