# ------------------------------------------------------ #
#### Prepare Art 17 checklists with recorded presence ####
# Here we load the Natura 2000 checklists of recorded species
# and habitats, linking them with the actual area of the sites
# Field explanations:
# https://www.eea.europa.eu/data-and-maps/data/natura-14/natura-2000-tabular-data-12-tables

# Get the layers and number
lyrs <- sf::st_layers("temporary_data/N2k.gpkg")

# Get the N2k sites
n2k <- sf::st_read("temporary_data/N2k.gpkg",layer = "NATURA2000SITES", quiet = TRUE) |>
  dplyr::filter(is.na(MARINE_AREA_PERCENTAGE) | MARINE_AREA_PERCENTAGE<5)
# Calculate the area from the actual shapes for comparison
# n2k_areas <- sf::st_read("temporary_data/N2k.gpkg", layer = "NaturaSite_polygon", quiet = TRUE)
# n2k_areas$area_ha <- units::set_units(sf::st_area(n2k_areas), "ha")
# ar <- n2k_areas |> dplyr::select(SITECODE, area_ha) |> sf::st_drop_geometry() |>
#   dplyr::group_by(SITECODE) |> dplyr::reframe(area_ha = sum(area_ha))
# n2k <- n2k |> dplyr::left_join(ar, by = "SITECODE")
# cor(n2k$area_ha,n2k$AREAHA) # 0.9998
# rm(ar,n2k_areas)

# Get the species and habitat layer
n2k_species <- sf::st_read("temporary_data/N2k.gpkg",layer = "SPECIES",quiet = TRUE)
n2k_habitats <- sf::st_read("temporary_data/N2k.gpkg",layer = "HABITATS",quiet = TRUE)

# Get list of annexes
check_species <- read.csv("temporary_data/Article17_2020_dataset_csv/Article17_2020_species_check_list.csv") |>
  dplyr::filter(group != "Fish")
n2k_species <- n2k_species |> dplyr::filter(SPECIESCODE %in% check_species$assessment_speciescode)
check_habitats <- read.csv("temporary_data/Article17_2020_dataset_csv/Article17_2020_habitats_check_list.csv")
n2k_habitats <- n2k_habitats |> dplyr::filter(HABITATCODE %in% check_habitats$habitatcode)

# --- #
# Remove non_presence site
# In case that a species no longer exists in the site a value "1" is entered.
n2k_species <- n2k_species |> dplyr::filter(NONPRESENCEINSITE %in% c(0,NA))
n2k_habitats <- n2k_habitats |> dplyr::filter(NON_PRESENCE_IN_SITE %in% c(0,NA))

# Filter to relevant columns
n2k_species <- n2k_species |>
  dplyr::select(COUNTRY_CODE, SITECODE,SPECIESNAME,SPECIESCODE,
                SPGROUP, LOWERBOUND, UPPERBOUND) |> distinct()

n2k_habitats <- n2k_habitats |>
  dplyr::select(COUNTRY_CODE, SITECODE,HABITATCODE) |>
  distinct()

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

# Get relevant columns from N2k
n2k <- n2k |> dplyr::select(COUNTRY_CODE,SITECODE,SITETYPE,AREAHA,DATE) |> distinct()

# For Area (in ha)
# The value -99 is given for sites for witch the areas is unknown.
# A value of 0 cab be correct if the site is a cave or cliff
#   (--> Enter 1 ha instead assuming a full 100m grid cell is covered)
n2k$AREAHA[n2k$AREAHA==0] <- 0.01
assertthat::assert_that(all(n2k$AREAHA>0))

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

#### Art 17 - Full species and distribution areas ####
# --- #
# Get area of ranges
lyr <- sf::st_layers("temporary_data/Art17_2018.gpkg")

# https://www.eea.europa.eu/data-and-maps/data/article-17-database-habitats-directive-92-43-eec-2/article-17-2015-dataset
distr_species <- sf::st_read("temporary_data/Art17_2018.gpkg",
                             layer = "Art17_species_distribution_2013_2018_EU",
                             quiet = TRUE)
distr_species <- distr_species |>
  dplyr::filter(maptype == "Distribution") |>
  dplyr::select(code, region) |>
  dplyr::distinct()

# Add new area in
distr_species$area_ha <- units::set_units(sf::st_area(distr_species),"ha")

# Summarize per code
ar_species <- distr_species |>
  sf::st_drop_geometry() |>
  dplyr::group_by(code) |>
  dplyr::reframe(distribution_ha = sum(area_ha,na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(distribution_ha = units::drop_units(distribution_ha))

# --- #
distr_habitats <- sf::st_read("temporary_data/Art17_2018.gpkg",
                              layer = "Art17_habitats_distribution_2013_2018_EU",
                              quiet = TRUE)
distr_habitats <- distr_habitats |>
  dplyr::filter(maptype == "Distribution") |>
  dplyr::select(habitatcode, region) |>
  dplyr::distinct()

# Add new area ina
distr_habitats$area_ha <- units::set_units(sf::st_area(distr_habitats),"ha")

# Summarize per code
ar_habitats <- distr_habitats |>
  sf::st_drop_geometry() |>
  dplyr::group_by(habitatcode) |>
  dplyr::reframe(distribution_ha = sum(area_ha)) |>
  dplyr::ungroup() |>
  dplyr::mutate(distribution_ha = units::drop_units(distribution_ha))

# ----------- #
#### Merge distribution and covered area ####
# Here we merge the N2k checklists with the full distribution data,
# also calculating representation as proportion covered within.

# First we merge the n2k data with the checklists
full1 <- n2k |> dplyr::left_join(n2k_habitats, by = c("COUNTRY_CODE", "SITECODE"))
full2 <- n2k |> dplyr::left_join(n2k_species, by = c("COUNTRY_CODE", "SITECODE"))

# Next we merge in the full distribution data by type
# This will be one to many relationship
full1 <- full1 |> dplyr::left_join(ar_habitats, by = c("HABITATCODE" = "habitatcode"))
full2 <- full2 |> dplyr::left_join(ar_species, by = c("SPECIESCODE" = "code"))

# Remove all for which we don't have full distribution matches
full1 <- full1 |> tidyr::drop_na(distribution_ha)
full2 <- full2 |> tidyr::drop_na(distribution_ha)

# What remains?
message("Habitats: ", n_distinct(full1$HABITATCODE), " in ", n_distinct(full1$SITECODE), " N2k sites")
message("Species: ", n_distinct(full2$SPECIESCODE), " in ", n_distinct(full2$SITECODE), " N2k sites")

# Add proportion
full1 <- full1 |> dplyr::rowwise() |>
  # Min is necessary as there 5 island sites with greater N2k areas than distributions
  dplyr::mutate(PROPCOVERED = min(AREAHA / distribution_ha,1))

full2 <- full2 |> dplyr::rowwise() |>
  # Min is necessary as there 5 sites with greater N2k areas than entire distributions
  dplyr::mutate(PROPCOVERED = min(AREAHA / distribution_ha,1))

# Correct those areas by setting their distribution to the AREAHA too
full1$distribution_ha <- ifelse(full1$AREAHA>full1$distribution_ha,full1$AREAHA,full1$distribution_ha)
full2$distribution_ha <- ifelse(full2$AREAHA>full2$distribution_ha,full2$AREAHA,full2$distribution_ha)

# Save the output as objects
saveRDS(full1, "temporary_data/n2k_habitatscovered.rds")
saveRDS(full2, "temporary_data/n2k_speciescovered.rds")

stop("Natura 2000 representation done and saved.")
