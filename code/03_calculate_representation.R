# # # # # # # # # # # # # # # #
# Purpose of this script is to calculate the feature representation by
# intersecting existing range data with protected area datasets.
#
# Author:
# Martin Jung (IIASA) - 2024
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
library(doParallel)
library(tidyr)
library(lubridate)

# Load custom functions
source("code/00_functions.R")

# Calculate from spatial data or from checklists?
spatial <- TRUE
aggregate_for_speed <- TRUE # Aggregate the original n2k to the Art17 grid. This increases the speed.

# TODO: Correct for land-use?
aoh <- FALSE

# Get list of conservation area layers
ll <- list.files("temporary_data/",full.names = TRUE)
ll <- ll[has_extension(ll,"tif")]
if(aggregate_for_speed) ll <- ll[grep("frac", ll)] else ll <- ll[grep("frac",ll,invert = TRUE)]
ll <- ll[grep("full", ll, invert = TRUE)]
consar <- terra::rast(ll)
names(consar) <- basename(tools::file_path_sans_ext(ll))
terra::time(consar) <- str_split(basename(tools::file_path_sans_ext(ll)),"_",simplify = TRUE)[,ifelse(aggregate_for_speed,3,2)] |> as.numeric()

if(spatial){
  # Load Biogeogregions for intersection with Art 12
  biogeo <- sf::st_read("temporary_data/BiogeoRegions2016.shp", quiet = TRUE) |>
    sf::st_transform(crs = sf::st_crs(consar)) |> dplyr::select(code,pre_2012) |>
    dplyr::rename("region" = "pre_2012") #|> sf::st_buffer(dist = 2000)
  # Rasterize to background
  bgbio <- terra::rasterize(biogeo, consar, field = "region",touches = TRUE)
  bgbio <- terra::focal(bgbio, w=15,na.policy = "only", fun = "modal")
  bgbio <- as.factor(bgbio)
  levels(bgbio) <- biogeo$region
  rm(biogeo)
}

# Checks
assertthat::assert_that(
  file.exists("temporary_data/N2k.gpkg"),
  file.exists("temporary_data/Art17_2018.gpkg")
)

if(!spatial) source("code/02a_calculate_representation_n2k.R")

# ------------------------------------------------------ #
#### Calculate total distribution and intersect with N2k data  ####
# Here we spatially intersect the conservation areas per
# time frame with the spatial distribution records of the species.
# We furthermore record the full distribution in the respective dataset

# Get sensitive species for subset
# NOTE: This is not a publicly available dataset!
path_art17 <- "temporary_data/Art12_Art17_SensitiveSpecies/eea_v_3035_10_mio_art17-2013-2018_i_2013-2018_v01_r00/Art17-2013-2018_GPKG/art17_2013_2018_internal.gpkg"
lyrs <- sf::st_layers(path_art17)

#registerDoParallel(cores = 10)

# Loop through each assessment type and save
for(name in lyrs$name[1:4]){ # name = lyrs$name[2]
  message(name)
  #if(file.exists(paste0("temporary_data/Speciesareas__",name, ".rds"))) next()

  data_sens <- sf::read_sf(path_art17,name,quiet = TRUE)

  # Output data.frame
  results <- data.frame()

  # Check that name is correct
  if(!utils::hasName(data_sens, "code")) data_sens$code <- data_sens$habitatcode
  if(!utils::hasName(data_sens, "code")) data_sens$code <- data_sens$speciescode
  if(!utils::hasName(data_sens, "country")) data_sens$country <- "All"

  assertthat::assert_that(utils::hasName(data_sens, "code"))

  pb <- progress::progress_bar$new(total = length(unique(data_sens$code)))
  # proc <- foreach(co = iter(data_sens$code),
  #         .packages = c("dplyr", "exactextractr","assertthat")) %dopar% {
  for(co in unique(data_sens$code)){ # co = unique(data_sens$code)[10]
    sub <- subset(data_sens, code == co)
    pb$tick()

    if(aoh) {
      # Apply to conservation areas by creating a subset and then applying
      # management specific priors using the reclassified Corine layer
      stop("TODO")
    }

    # For year in stack
    for(lyr in names(consar)){ # names(consar)[4]->lyr
      # Exact extract
      ex <- exactextractr::exact_extract(consar[[lyr]], sub, fun = "sum",
                                         force_df = TRUE, progress = FALSE)
      # Combine with sub and summarize
      o <- dplyr::bind_cols(sub, ex) |>
        dplyr::group_by(maptype, category, region, country, code) |>
        dplyr::reframe(conservedarea_m2 =
                         units::set_units(
                           sum(sum) * Reduce("*",terra::res(consar)), "m2"
                         )
        )
      o <- o |> dplyr::mutate(year = terra::time( consar[[lyr]] ), .before = "code")
      # And name
      o <- o |> dplyr::mutate(name = name, .before = "maptype")

      # Also add total distribution area
      o <- dplyr::left_join(
        o,
        sub |> dplyr::group_by(maptype, category, region, country, code) |>
          dplyr::reframe(totalarea_m2 = sum(sf::st_area(geom)))
      )
      # Check that all area equal/smaller
      assertthat::assert_that(
        all(o$totalarea_m2>=o$conservedarea_m2)
      )

      # return(o)
      # Add
      results <- bind_rows(results, o)
      rm(o)
    }

  }

  # Save the results
  saveRDS(object = results, file = paste0("temporary_data/Speciesareas__",name, ".rds"))
  rm(results)
}

# --- #
# Article 12
path_art12 <- "temporary_data/Art12_Art17_SensitiveSpecies/eea_v_3035_10_mio_art12-2013-2018_i_2013-2018_v01_r00/Art12-2013-2018_GPKG/ART12_3035_distribution_data_with_sensitive.gpkg"
lyrs <- sf::st_layers(path_art12)

# Loop through each assessment type and save
for(name in lyrs$name){ # name = lyrs$name
  message(name)
  # if(file.exists(paste0("temporary_data/Speciesareas__",name, ".rds"))) next()

  data_sens <- sf::read_sf(path_art12,name,quiet = TRUE)

  # Output data.frame
  results <- data.frame()

  # Check that name is correct
  if(!utils::hasName(data_sens, "code")) data_sens$code <- data_sens$speciescode
  if(!utils::hasName(data_sens, "country")) data_sens$country <- data_sens$country

  # Art 12 has no region, thus we intersect it to determine the region
  if(!utils::hasName(data_sens, "region")){
    data_sens$region <- exactextractr::exact_extract(bgbio, data_sens, fun = "mode", summarize_df = F)
    # Recodify based on the factor levels
    data_sens$region <- levels(bgbio)[[1]][,2][match(data_sens$region, levels(bgbio)[[1]][,1])]
    # For missing ones (coastal cell), take nearest region
    assertthat::assert_that(!anyNA(data_sens$region))
    rm(ex,check)
  }
  assertthat::assert_that(utils::hasName(data_sens, "code"))

  pb <- progress::progress_bar$new(total = length(unique(data_sens$code)))
  for(co in unique(data_sens$code)){
    sub <- subset(data_sens, code == co)
    pb$tick()

    if(aoh) {
      # Apply to conservation areas by creating a subset and then applying
      # management specific priors using the reclassified Corine layer
      stop("TODO")
    }

    # For year in stack
    for(lyr in names(consar)){
      # Exact extract
      ex <- exactextractr::exact_extract(consar[[lyr]], sub, fun = "sum",
                                         force_df = TRUE)
      # Combine with sub and summarize
      o <- dplyr::bind_cols(sub, ex) |>
        dplyr::group_by(maptype, season, region, country, code) |>
        dplyr::reframe(conservedarea_m2 =
                         units::set_units(
                           sum(sum) * Reduce("*",terra::res(consar)), "m2"
                         )
        )
      o <- o |> dplyr::mutate(year = terra::time( consar[[lyr]] ), .before = "code")
      # And name
      o <- o |> dplyr::mutate(name = name, .before = "maptype")
      # And other species_code
      o <- o |> dplyr::mutate(speciescodeEU = unique(sub$speciescodeEU), .before = "code")

      # Also add total distribution area
      o <- dplyr::left_join(
        o,
        sub |> dplyr::group_by(maptype, season, region, country, code) |>
          dplyr::reframe(totalarea_m2 = sum(sf::st_area(geom)))
      )
      # Check that all area equal/smaller
      assertthat::assert_that(
        all(o$totalarea_m2>=o$conservedarea_m2)
      )

      # Add
      results <- bind_rows(results, o)
      rm(o)
    }
  }

  # Save the results
  saveRDS(object = results, file = paste0("temporary_data/Speciesareas__",name, ".rds"))
  rm(results)
}

stop("Spatial: Natura 2000 representation done and saved.")
