# # # # # # # # # # # # # # # #
# Purpose of this script is to calculate the rasterize the protection layers
# to the spatial grain and extent of the Corine accounting layer.
#
# Author:
# Martin Jung (IIASA) - 2024
# # # # # # # # # # # # # # # #

library(terra)
library(sf)
library(ibis.iSDM)
library(assertthat)
# Load custom functions
source("code/00_functions.R")

# Direction to temporary folder
path_temp <- "temporary_data/"
dir.create(path_temp, showWarnings = FALSE)

# Path to Corine accounting layer for 2000
path_corine <- paste0(path_temp, "CLC_account/")
assertthat::assert_that(dir.exists(path_corine))

# Path to protected area layers
path_n2k <- paste0(path_temp, "/N2k.gpkg")
assertthat::assert_that(file.exists(path_n2k))

# Path NUTS
path_nuts <- paste0(path_temp, "NUTS_RG_20M_2024_4326.shp/NUTS_RG_20M_2024_4326.shp")
assertthat::assert_that(file.exists(path_nuts))

# ----------------- #
#### Create background and rasterize conservation areas ####
# Here we create a binary background that contains the conservation areas
# subset used for the accounting of representation.

corine <- terra::rast( paste0(path_corine, "CLC2000ACC_V2018_20.tif") )
background <- emptyraster(corine)

# Also get the reference grid
background_10km <- terra::rast("../11.1.2_ReferenceGrids/ReferenceGrid_Europe_bin_10000m.tif") |>
  terra::project(terra::crs(background))

# Get the N2k areas
lyrs <- sf::st_layers("temporary_data/N2k.gpkg")

# Get the N2k sites
n2k <- sf::st_read("temporary_data/N2k.gpkg",layer = "NATURA2000SITES", quiet = TRUE) |>
  dplyr::filter(is.na(MARINE_AREA_PERCENTAGE) | MARINE_AREA_PERCENTAGE<5)

# Get the layers and number
n2k_sites <- sf::st_read("temporary_data/N2k.gpkg",layer = "NaturaSite_polygon", quiet = TRUE)
n2k_sites <- n2k_sites |> dplyr::filter(SITECODE %in% n2k$SITECODE)

# --- #
# Classify years
# Determine new date dependent on type and availability
# A: SPAs (Special Protection Areas - sites designated under the Birds Directive);
# B: SCIs and SACs (Sites of Community Importance and Special Areas of Conservation - sites designated under the Habitats Directive);
# C: where SPAs and SCIs/SACs boundaries are identical (sites designated under both directives).
n2k$DATE <- ifelse(is.na(n2k$DATE_SPA),
                   ifelse(is.na(n2k$DATE_SAC),
                          ifelse(is.na(n2k$DATE_PROP_SCI),
                                 # Assume the compilation date as the earliest
                                 # Date of recognition
                                 as.character(n2k$DATE_COMPILATION),
                                 as.character(n2k$DATE_PROP_SCI)),
                          as.character(n2k$DATE_SAC)),
                   as.character(n2k$DATE_SPA))

# Get missing dates (France only) and assign DATE_CONF_SCI here
n2k[which(is.na(n2k$DATE)),"DATE"] <- as.character(n2k[which(is.na(n2k$DATE)),"DATE_CONF_SCI"])

# For all remaining ones assign earliest date as date of N2k enactment
n2k[which(is.na(n2k$DATE)),"DATE"] <- as.character("1992-05-01")
assertthat::assert_that(!anyNA(n2k$DATE))

# Format and categorize the dates
n2k$DATE <- as.Date(n2k$DATE)
assertthat::assert_that(!anyNA(n2k$DATE))
# Correct that one wrong entry in Bulgaria (typo to 2013)
n2k[which(n2k$DATE=="2031-03-01"),"DATE"] <- as.Date("2013-03-01")
assertthat::assert_that(!anyNA(n2k$DATE))

# Add a year and group
n2k$YEAR <- lubridate::year(n2k$DATE)

n2k$YEAR_GROUP <- cut(n2k$YEAR, c(1900,2000,2006,2012,2018,2024),
                      labels = c("<2000", "2000-2006", "2006-2012",
                                 "2012-2018", "2018-2024"),include.lowest = TRUE)
assertthat::assert_that(!anyNA(n2k$YEAR_GROUP))

# --- #
#### Rasterize for each year group ####
# Here the idea is write one layer per YEAR GROUP
# so as to capture expansion per years. Note

for(yg in levels(n2k$YEAR_GROUP)){ # levels(n2k$YEAR_GROUP)[2]->yg
  message(yg)
  # Subset to all year prior to the last as cumulative expansion
  lastyear <- switch (as.character(yg),
    "<2000" = 2000,
    "2000-2006" = 2006,
    "2006-2012" = 2012,
    "2012-2018" = 2018,
    "2018-2024" = 2024
  )
  if(file.exists(paste0(path_temp, "natura2000_",as.character(lastyear),".tif"))) next()

  # Subset including all sites younger and including this year
  sub <- subset(n2k, YEAR <= lastyear)

  n2k_sub <- n2k_sites |> dplyr::filter(SITECODE %in% unique(sub$SITECODE))
  # Rasterize and write the whole output
  n2k_full <- terra::rasterize(n2k_sub, background, field = 1, touches = TRUE)

  terra::writeRaster(n2k_full,
                     filename = paste0(path_temp, "natura2000_",as.character(lastyear),".tif"),
                     datatype = 'INT1U',
                     overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE"),
                     progress = TRUE
  )

  # Also for 10k
  n2k_full <- terra::rasterize(n2k_sub, background_10km, cover = TRUE, touches = TRUE)
  terra::writeRaster(n2k_full,
                     filename = paste0(path_temp, "natura2000_frac_",as.character(lastyear),".tif"),
                     datatype = 'FLT4S',
                     overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE"),
                     progress = TRUE
  )
  rm(n2k_full);gc()
}

# ------------------- #
#### Calculate area expansion by time ####
# Here we simply calculate the total area by timeslot
ll <- list.files("temporary_data/",full.names = TRUE)
ll <- ll[has_extension(ll,"tif")]
ll <- ll[grep("frac", ll)] #else ll <- ll[grep("frac",ll,invert = TRUE)]
ll <- ll[grep("full", ll, invert = TRUE)]
consar <- terra::rast(ll)
names(consar) <- basename(tools::file_path_sans_ext(ll))
terra::time(consar) <- str_split(basename(tools::file_path_sans_ext(ll)),"_",simplify = TRUE)[,3] |> as.numeric()

# Also load the NUTS file
nuts <- sf::st_read(path_nuts, quiet = TRUE) |> dplyr::filter(LEVL_CODE == 0) |>
  # Exclude non-covered countries
  dplyr::filter(!(CNTR_CODE %in% c("TR","NO","CH", "IS", "UA","XK"))) |>
  # Add area
  dplyr::mutate(area_km2 = units::set_units(sf::st_area(geometry), "km2"))
#plot(nuts['NUTS_ID'])

# Aggregate overall
n2k <- terra::global(consar * terra::cellSize(consar[[1]], unit = "km"),
                     'sum', na.rm = TRUE) |> tibble::rownames_to_column("year") |>
  dplyr::rename(area_km2 = sum) |>
  dplyr::mutate(year = stringr::str_split(year, "_",simplify = TRUE)[,3],
                area_km2 = units::set_units(area_km2,"km2"))
n2k$year <- factor(n2k$year,
                  levels = c(2000, 2006, 2012, 2018, 2024),
                  labels = c("<2000", "2000-2006", "2006-2012","2012-2018", "2018-2024"))
# Finally add overall proportion
n2k$fullprop <- n2k$area_km2 / sum(nuts$area_km2)

# Save output
write.csv(n2k, "export/Overall_n2k.csv",row.names = FALSE)
