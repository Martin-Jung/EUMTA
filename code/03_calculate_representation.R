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
library(ibis.iSDM)

# Load custom functions
source("code/00_functions.R")

# Calculate from spatial data or from checklists?
spatial <- TRUE
aggregate_for_speed <- TRUE # Aggregate the original n2k to the Art17 grid. This increases the speed.

cores <- 10

# TODO: Correct for land-use?
aoh <- FALSE

# --- #

# Get list of conservation area layers
ll <- list.files("temporary_data/",full.names = TRUE)
ll <- ll[has_extension(ll,"tif")]
if(aggregate_for_speed) ll <- ll[grep("frac", ll)] else ll <- ll[grep("frac",ll,invert = TRUE)]
ll <- ll[grep("full", ll, invert = TRUE)]

tt <- stringr::str_split(basename(ll),"_",simplify = TRUE)[,1]
# Now make a list per type
cas <- list()
for(ctype in unique(tt)){
  consar <- terra::rast(ll[which(tt==ctype)])
  names(consar) <- basename(tools::file_path_sans_ext(ll[which(tt==ctype)]))
  terra::time(consar) <- str_split(basename(tools::file_path_sans_ext(ll[which(tt==ctype)])),"_",simplify = TRUE)[,ifelse(aggregate_for_speed,3,2)] |> as.numeric()
  cas[[ctype]] <- consar
  rm(consar)
}
cas <- Reduce("c", cas) # Combine all

if(spatial){
  # Load Biogeogregions for intersection with Art 12
  biogeo <- sf::st_read("temporary_data/BiogeoRegions2016.shp", quiet = TRUE) |>
    sf::st_transform(crs = sf::st_crs(cas[[1]])) |> dplyr::select(code,pre_2012) |>
    dplyr::rename("region" = "pre_2012") #|> sf::st_buffer(dist = 2000)
  # Rasterize to background
  bgbio <- terra::rasterize(biogeo, cas[[1]], field = "region",touches = TRUE)
  bgbio <- terra::focal(bgbio, w=15,na.policy = "only", fun = "modal")
  bgbio <- as.factor(bgbio)
  levels(bgbio) <- biogeo$region
  rm(biogeo)
}

# Multiply cas with cellsize in km2
cas <- cas * terra::cellSize(cas, unit="km")

# Checks
assertthat::assert_that(
  is.numeric(cores),
  is.Raster(cas), terra::nlyr(cas)>4,
  file.exists("temporary_data/N2k.gpkg"),
  file.exists("temporary_data/Art17_2018.gpkg")
)

if(!spatial) source("code/02a_calculate_representation_n2k.R")

# ------------------------------------------------------ #
#### Calculate total distribution and intersect with rasterized conservation data  ####
# Here we spatially intersect the conservation areas per
# time frame with the spatial distribution records of the species.
# We furthermore record the full distribution in the respective dataset

# Get sensitive species for subset
# NOTE: This is not a publicly available dataset!
path_art17 <- "temporary_data/Art12_Art17_SensitiveSpecies/eea_v_3035_10_mio_art17-2013-2018_i_2013-2018_v01_r00/Art17-2013-2018_GPKG/art17_2013_2018_internal.gpkg"
lyrs <- sf::st_layers(path_art17)

# Set parallel processing
# cl <- parallel::makeCluster(cores)
# registerDoParallel(cl = cl)

# Loop through each assessment type focussing on EU
for(name in lyrs$name[1:4]){ # name = lyrs$name[3]
  message(name)
  if(file.exists(paste0("temporary_data/Speciesareas__",name, ".rds"))) next()

  data_sens <- sf::read_sf(path_art17,name,quiet = TRUE)

  # Check that name is correct
  if(!utils::hasName(data_sens, "code")) data_sens$code <- data_sens$habitatcode
  if(!utils::hasName(data_sens, "code")) data_sens$code <- data_sens$speciescode
  if(!utils::hasName(data_sens, "country")) data_sens$country <- "All"

  assertthat::assert_that(utils::hasName(data_sens, "code"))

  results <- foreach(co = unique(data_sens$code),
                  .inorder = FALSE,
                  .combine = rbind,
                  .multicombine = TRUE,
                  .export = c("data_sens", "cas", "aoh"),
                  .packages = c("dplyr", "exactextractr","terra","assertthat")
                  ) %do% {
    # pb <- progress::progress_bar$new(total = length(unique(data_sens$code)))
    # for(co in unique(data_sens$code)){ # co = unique(data_sens$code)[1]
    sub <- subset(data_sens, code == co)
    # pb$tick()

    if(aoh) {
      # Apply to conservation areas by creating a subset and then applying
      # management specific priors using the reclassified Corine layer
      stop("TODO")
    }

    # Exact extract
    ex <- exactextractr::exact_extract(cas,
                                       sub, fun = "sum",
                                       force_df = TRUE, progress = FALSE)
    names(ex) <- stringr::str_remove(names(ex), "sum.")
    assertthat::assert_that(nrow(sub) == nrow(ex))

    # Combine with sub and summarize
    o <- dplyr::bind_cols(sub, ex)
    # Split and clear value
    o <- o |> pivot_longer(cols= names(ex)) |>
      tidyr::separate(name, into = c("ctype","sep","year"),sep = "_") |>
      sf::st_drop_geometry() |>
      dplyr::group_by(maptype, category, ctype, region,country, year, code) |>
      dplyr::reframe(conservedarea_km2 = units::set_units(sum(value),"km2") )
    # Add dataset name
    o <- o |> dplyr::mutate(name = name, .before = "maptype")

    # Also add total distribution area
    o <- dplyr::left_join(
      o,
      sub |> dplyr::group_by(maptype, category, region, country, code) |>
        dplyr::reframe(totalarea_km2 = units::set_units( sum(sf::st_area(geom)),"km2") )
    )

    # Check for rounding issues with area distortion in which case, take the maximum
    if(any(o$conservedarea_km2>=o$totalarea_km2)){
      message("Area rounding issue found for ", co)
      o <- o |> dplyr::rowwise() |>
        dplyr::mutate(totalarea_km2 = pmax(totalarea_km2,conservedarea_km2))
    }
    # Check that all area equal/smaller
    assertthat::assert_that(
      is.numeric(o$conservedarea_km2),is.numeric(o$totalarea_km2),
      all(o$totalarea_km2>=o$conservedarea_km2)
    )
    return(o)
    rm(o)
  }

  assertthat::assert_that(is.data.frame(results))
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
  if(file.exists(paste0("temporary_data/Speciesareas__",name, ".rds"))) next()

  data_sens <- sf::read_sf(path_art12,name,quiet = TRUE)

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

  results <- foreach(co = unique(data_sens$code),
                     .inorder = FALSE,
                     .combine = rbind,
                     .multicombine = TRUE,
                     .export = c("data_sens", "cas", "aoh"),
                     .packages = c("dplyr", "exactextractr","terra","assertthat")
  ) %do% {
    # pb <- progress::progress_bar$new(total = length(unique(data_sens$code)))
    # for(co in unique(data_sens$code)){ # co = unique(data_sens$code)[10]
    sub <- subset(data_sens, code == co)
    # pb$tick()

    if(aoh) {
      # Apply to conservation areas by creating a subset and then applying
      # management specific priors using the reclassified Corine layer
      stop("TODO")
    }

    # Exact extract
    ex <- exactextractr::exact_extract(cas,
                                       sub, fun = "sum",
                                       force_df = TRUE, progress = FALSE)
    names(ex) <- stringr::str_remove(names(ex), "sum.")
    # Combine with sub and summarize
    o <- dplyr::bind_cols(sub, ex)
    # Split and clear value
    o <- o |> pivot_longer(cols= names(ex)) |>
      tidyr::separate(name, into = c("ctype","sep","year"),sep = "_") |>
      dplyr::group_by(maptype, season, ctype, region, country, year, code) |>
      dplyr::reframe(conservedarea_km2 = units::set_units(sum(value),"km2") )
    # Add dataset name
    o <- o |> dplyr::mutate(name = name, .before = "maptype")
    # And other species_code
    o <- o |> dplyr::mutate(speciescodeEU = unique(sub$speciescodeEU), .before = "code")

    # Also add total distribution area
    o <- dplyr::left_join(
      o,
      sub |> dplyr::group_by(maptype, season, region, country, code) |>
        dplyr::reframe(totalarea_km2 = units::set_units( sum(sf::st_area(geom)),"km2") )
    )
    # Check for rounding issues with area distortion in which case, take the maximum
    if(any(o$conservedarea_km2>=o$totalarea_km2)){
      message("Area rounding issue found for ", co)
      o <- o |> dplyr::rowwise() |>
        dplyr::mutate(totalarea_km2 = pmax(totalarea_km2,conservedarea_km2))
    }
    # Check that all area equal/smaller
    assertthat::assert_that(
      is.numeric(o$conservedarea_km2),is.numeric(o$totalarea_km2),
      all(o$totalarea_km2>=o$conservedarea_km2)
    )
    return(o)
    rm(o)
  }

  assertthat::assert_that(is.data.frame(results))

  # Save the results
  saveRDS(object = results, file = paste0("temporary_data/Speciesareas__",name, ".rds"))
  rm(results)
}

stop("Spatial: Natura 2000 representation done and saved.")
