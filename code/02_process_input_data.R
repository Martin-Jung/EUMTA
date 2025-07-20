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
library(lubridate)
library(dplyr)
library(stringr)
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

# Path to CDDA layer
path_cdda <- paste0(path_temp, "/CDDA.gpkg")
assertthat::assert_that(file.exists(path_cdda))

# Path NUTS
path_nuts <- paste0(path_temp, "NUTS_RG_20M_2024_4326.shp/NUTS_RG_20M_2024_4326.shp")
assertthat::assert_that(file.exists(path_nuts))

# Create background and rasterize conservation areas #
corine <- terra::rast( paste0(path_corine, "CLC2000ACC_V2018_20.tif") )
background <- emptyraster(corine)

# Also get the reference grid
background_10km <- terra::rast("../11.1.2_ReferenceGrids/ReferenceGrid_Europe_bin_10000m.tif") |>
  terra::project(terra::crs(background))

# ----------------- #
#### Process Natura 2000 data ####
# Here we create a binary background that contains the conservation areas
# subset used for the accounting of representation.

# Get the N2k areas
lyrs <- sf::st_layers(path_n2k)

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

## Old approach
# n2k$DATE <- ifelse(is.na(n2k$DATE_SPA),
#                    ifelse(is.na(n2k$DATE_SAC),
#                           ifelse(is.na(n2k$DATE_PROP_SCI),
#                                  # Assume the compilation date as the earliest
#                                  # Date of recognition
#                                  as.character(n2k$DATE_COMPILATION),
#                                  as.character(n2k$DATE_PROP_SCI)),
#                           as.character(n2k$DATE_SAC)),
#                    as.character(n2k$DATE_SPA))

# New query
n2k$A_year <- ifelse(n2k$SITETYPE == 'A', lubridate::year(n2k$DATE_SPA), 0)
n2k$B_year <- ifelse(n2k$SITETYPE == 'B', lubridate::year(n2k$DATE_PROP_SCI), 0)
n2k$C_year <- ifelse(n2k$SITETYPE == 'C',
                     ifelse(lubridate::year(n2k$DATE_PROP_SCI) < lubridate::year(n2k$DATE_SPA),
                            lubridate::year(n2k$DATE_PROP_SCI), lubridate::year(n2k$DATE_SPA)),
                     0
)

# Assign the date as the non-zero year
n2k$DATE <- ifelse(n2k$A_year > 0, n2k$A_year,
                        ifelse(n2k$B_year > 0, n2k$B_year,
                               ifelse(n2k$C_year > 0, n2k$C_year, NA)))

# How many have an existing date?
message("Number of N2k sites with a date: ",
        nrow(n2k[which(!is.na(n2k$DATE)),]),
        " out of ", nrow(n2k), " (",
        round(nrow(n2k[which(!is.na(n2k$DATE)),]) / nrow(n2k) * 100, 1), "%)")

table(n2k$COUNTRY_CODE[which(is.na(n2k$DATE))]) # Almost all french N2k sites here...

# Get missing dates (almost exclusively France only) and assign DATE_CONF_SCI here
n2k[which(is.na(n2k$DATE)),"DATE"] <- lubridate::year( as.character(n2k[which(is.na(n2k$DATE)),"DATE_CONF_SCI"]) )

# For all remaining ones assign earliest date as date of N2k enactment
n2k[which(is.na(n2k$DATE)),"DATE"] <- lubridate::year(as.character("1992-05-01"))
assertthat::assert_that(!anyNA(n2k$DATE))

# Add a year and group
n2k$YEAR <- (n2k$DATE)

n2k$YEAR_GROUP <- cut(n2k$YEAR, c(1900,2000,2006,2012,2018,2024),
                      labels = c("<2000", "2000-2006", "2006-2012",
                                 "2012-2018", "2018-2024"),include.lowest = TRUE)
assertthat::assert_that(!anyNA(n2k$YEAR_GROUP))

#### Process CDDA data ####
# Here we create a binary background that contains the national designated areas
# subset used for the accounting of representation.

# Get the CDDA areas
lyrs <- sf::st_layers(path_cdda)

# Get the CDDA sites
cdda <- sf::st_read(path_cdda,layer = "DesignatedArea", quiet = TRUE) |>
  # Remove marine areas, but retain border regions such as marineAndTerrestrial
  dplyr::filter(majorEcosystemType != "marine") |>
  dplyr::filter(majorEcosystemType == "terrestrial" | is.na(marineAreaPercentage))

# Get the layers and remove the ones not selected above
cdda_sites <- sf::st_read(path_cdda,layer = "ProtectedSite", quiet = TRUE)
cdda_sites <- cdda_sites |> dplyr::filter(cddaId %in% cdda$cddaId)
sf::st_geometry(cdda_sites) <- "geometry"
# Cast to multipolygon and to background
cdda_sites <- cdda_sites |> sf::st_transform(crs = sf::st_crs(background))

# --- Ignore point data as this is not compatable with the reference grid ---
cdda_sites <- cdda_sites |> dplyr::filter(geomType == "polygon")
# --- #
cdda_sites <- cdda_sites |> sf::st_cast('MULTIPOLYGON') # Cast to multypolygon

message("Number of cdda sites: ", dplyr::n_distinct(cdda_sites$cddaId))

# --- #
# Classify years
cdda$YEAR <- cdda$legalFoundationDate
# For all remaining ones assign earliest date as date of N2k enactment
cdda[which(is.na(cdda$YEAR)),"YEAR"] <- as.character("1992")
cdda$YEAR <- as.numeric(cdda$YEAR)
assertthat::assert_that(!anyNA(cdda$YEAR),
                        all(cdda$YEAR <= lubridate::year(lubridate::today()) ))

cdda$YEAR_GROUP <- cut(cdda$YEAR, c(min(cdda$YEAR),2000,2006,2012,2018,2024),
                      labels = c("<2000", "2000-2006", "2006-2012",
                                 "2012-2018", "2018-2024"),include.lowest = TRUE)
assertthat::assert_that(!anyNA(cdda$YEAR_GROUP))

# --- #
#### Rasterize for each year group ####
# Here the idea is write one layer per YEAR GROUP
# so as to capture expansion per years across the reference grid.

# First for Natura2000
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

# Then for CDDA
for(yg in levels(cdda$YEAR_GROUP)){ # levels(cdda$YEAR_GROUP)[2]->yg
  message(yg)
  # Subset to all year prior to the last as cumulative expansion
  lastyear <- switch (as.character(yg),
                      "<2000" = 2000,
                      "2000-2006" = 2006,
                      "2006-2012" = 2012,
                      "2012-2018" = 2018,
                      "2018-2024" = 2024
  )
  if(file.exists(paste0(path_temp, "cdda_",as.character(lastyear),".tif"))) next()

  # Subset including all sites younger and including this year
  sub <- subset(cdda, YEAR <= lastyear)

  cdda_sub <- cdda_sites |> dplyr::filter(cddaId %in% unique(sub$cddaId))
  # Rasterize and write the whole output for the background
  cdda_full <- terra::rasterize(cdda_sub, background, field = 1, touches = TRUE)

  terra::writeRaster(cdda_full,
                     filename = paste0(path_temp, "cdda_",as.character(lastyear),".tif"),
                     datatype = 'INT1U',
                     overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE"),
                     progress = TRUE
  )

  # Also for 10k
  cdda_full <- terra::rasterize(cdda_sub, background_10km, cover = TRUE, touches = TRUE)
  terra::writeRaster(cdda_full,
                     filename = paste0(path_temp, "cdda_frac_",as.character(lastyear),".tif"),
                     datatype = 'FLT4S',
                     overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE"),
                     progress = TRUE
  )
  rm(cdda_full,n2k_full);gc()
}

# Then for both combined
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
  if(file.exists(paste0(path_temp, "combined_",as.character(lastyear),".tif")) ) next()

  # Subset including all sites younger and including this year
  sub1 <- subset(n2k, YEAR <= lastyear)
  sub2 <- subset(cdda, YEAR <= lastyear)

  # Filter
  n2k_sub <- n2k_sites |> dplyr::filter(SITECODE %in% unique(sub1$SITECODE)) |>
    sf::st_transform(crs = sf::st_crs(background))
  sf::st_geometry(n2k_sub) <- "geometry"
  cdda_sub <- cdda_sites |> dplyr::filter(cddaId %in% unique(sub2$cddaId)) |>
    sf::st_transform(crs = sf::st_crs(background))
  sf::st_geometry(cdda_sub) <- "geometry"

  # Then select relevant columns and rename and combine
  comb <- bind_rows(
    n2k_sub |> dplyr::rename(siteid = SITECODE) |> dplyr::select(siteid) |>
      dplyr::mutate(siteid = as.character(siteid)),
    cdda_sub |> dplyr::rename(siteid = cddaId) |> dplyr::select(siteid) |>
      dplyr::mutate(siteid = as.character(siteid))
  )

  # Rasterize and write the whole output
  comb_full <- terra::rasterize(comb, background, field = 1, touches = TRUE)

  terra::writeRaster(comb_full,
                     filename = paste0(path_temp, "combined_",as.character(lastyear),".tif"),
                     datatype = 'INT1U',
                     overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE"),
                     progress = TRUE
  )

  # Also for 10k
  comb_full <- terra::rasterize(comb, background_10km, cover = TRUE, touches = TRUE)
  terra::writeRaster(comb_full,
                     filename = paste0(path_temp, "combined_frac_",as.character(lastyear),".tif"),
                     datatype = 'FLT4S',
                     overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE"),
                     progress = TRUE
  )
  rm(comb_full, comb);gc()
}

message("Rasterized layers to reference grids!")

#### Calculate area expansion statistics by time ####
# Load the NUTS file
nuts <- sf::st_read(path_nuts, quiet = TRUE) |> dplyr::filter(LEVL_CODE == 0) |>
  # Exclude non-covered countries
  dplyr::filter(!(CNTR_CODE %in% c("TR","NO","CH", "IS", "UA","XK"))) |>
  # Add area
  dplyr::mutate(area_km2 = units::set_units(sf::st_area(geometry), "km2"))
#plot(nuts['NUTS_ID'])

# Here we simply calculate the total area by timeslot
ll <- list.files("temporary_data/",full.names = TRUE)
ll <- ll[has_extension(ll,"tif")]
ll <- ll[grep("frac", ll)] #else ll <- ll[grep("frac",ll,invert = TRUE)]
ll <- ll[grep("full", ll, invert = TRUE)]

# Now iterate through types
out <- data.frame()
tt <- stringr::str_split(basename(ll),"_",simplify = TRUE)[,1]
for(ctype in unique(tt)){
  message(ctype)

  consar <- terra::rast(ll[which(tt==ctype)])
  names(consar) <- basename(tools::file_path_sans_ext(ll[which(tt==ctype)]))
  terra::time(consar) <- str_split(basename(tools::file_path_sans_ext(ll[which(tt==ctype)])),"_",simplify = TRUE)[,3] |> as.numeric()

  # Rasterize the nuts area and then summarize
  ras_nuts <- terra::rasterize(nuts |> sf::st_transform(crs = sf::st_crs(consar)),
                               consar, field = 1)

  # Mask consar to be sure to cover only EU27
  consar <- terra::mask(consar, ras_nuts)

  # Aggregate overall
  ca <- terra::global(consar * terra::cellSize(consar[[1]], unit = "km"),
                       'sum', na.rm = TRUE) |> tibble::rownames_to_column("year") |>
    dplyr::rename(area_km2 = sum) |>
    dplyr::mutate(ctype = ctype,
                  year = stringr::str_split(year, "_",simplify = TRUE)[,3],
                  area_km2 = units::set_units(area_km2,"km2"))
  ca$year <- factor(ca$year,
                     levels = c(2000, 2006, 2012, 2018, 2024),
                     labels = c("<2000", "2000-2006", "2006-2012","2012-2018", "2018-2024"))

  # Finally add overall proportion
  ca$fullprop <- ca$area_km2 / terra::global(terra::cellSize(ras_nuts,unit="km"), "sum")[,1]*10

  out <- dplyr::bind_rows(out, ca)
  rm(ca,consar)
}

# Save output
write.csv(out, "export/Overall_areastatistics.csv",row.names = FALSE)
