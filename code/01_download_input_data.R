# # # # # # # # # # # # # # # #
# Purpose of this script is to download or respectively process the
# original spatial layers from the EEA. We will store those files in a temporary
# data folder as they can be quite large.
#
# Author:
# Martin Jung (IIASA) - 2024
# # # # # # # # # # # # # # # #

# ---- #
# Packages
library(tidyverse)
library(assertthat)
library(curl)

# Direction to temporary folder
path_temp <- "temporary_data/"
dir.create(path_temp, showWarnings = FALSE)

# Biogeographic regions
# Source: https://sdi.eea.europa.eu/catalogue/srv/eng/catalog.search#/metadata/c6d27566-e699-4d58-a132-bbe3fe01491b
url <- "https://sdi.eea.europa.eu/webdav/datastore/public/eea_v_3035_1_mio_biogeo-regions_p_2016_v01_r00/BiogeoRegions2016.shp"
curl::curl_download(url, destfile = paste0(path_temp, "BiogeoRegions.shp"))

# Get CDDA geopackage file from EEA for Europe #
# Source: https://sdi.eea.europa.eu/catalogue/srv/eng/catalog.search#/metadata/616ef48f-7196-4e30-b201-6c97808fa68a
url <- "https://sdi.eea.europa.eu/webdav/datastore/public/eea_v_3035_100_k_natda-poly_p_2023-2024_v22_r00/GPKG/NatDA_2024_v01_public_EuropeEPSG3035.gpkg"
curl::curl_download(url, destfile = paste0(path_temp, "CDDA.gpkg"))

# Get Natura 2000 data #
# Source: https://sdi.eea.europa.eu/catalogue/srv/eng/catalog.search#/metadata/95e717d4-81dc-415d-a8f0-fecdf7e686b0
url <- "https://sdi.eea.europa.eu/webdav/datastore/public/eea_v_3035_100_k_natura2000_p_2022_v01_r00/Natura2000_end2022.gpkg"
curl::curl_download(url, destfile = paste0(path_temp, "N2k.gpkg"))

# --- #
## Species and habitats ##

# Get Article 17 (Public!) - 2018
# Source: https://sdi.eea.europa.eu/catalogue/srv/eng/catalog.search#/metadata/9f71b3e3-f8ec-442b-a2d5-c3c190605ac4
url <- "https://sdi.eea.europa.eu/webdav/datastore/public/eea_v_3035_10_mio_art17-2013-2018_p_2013-2018_v01_r00/Art17-2013-2018_GPKG/art17_2013_2018_public.gpkg"
curl::curl_download(url, destfile = paste0(path_temp, "Art17_2018.gpkg"))

# Get Article 17 (Public!) - 2012
# Source: https://sdi.eea.europa.eu/catalogue/srv/eng/catalog.search#/metadata/7fc458f8-40e1-4528-87bf-62ef0896fbb3
url <- "https://sdi.eea.europa.eu/webdav/datastore/public/eea_v_3035_10_mio_art17-2007-2012_p_2007-2012_v01_r01/Art17-2007-2012_GPKG/ART17_2007_2012_public.gpkg"
curl::curl_download(url, destfile = paste0(path_temp, "Art17_2012.gpkg"))

# Get Article 12 (Public!) - 2018
# Source: https://sdi.eea.europa.eu/catalogue/srv/eng/catalog.search#/metadata/e2face16-f352-4aff-9e4f-0ad1306f89b5
url <- "https://sdi.eea.europa.eu/webdav/datastore/public/eea_v_3035_10_mio_art12-2013-2018_p_2013-2018_v01_r01/Art12-2013-2018_GPKG/ART12_3035_distribution_data_without_sensitive.gpkg"
curl::curl_download(url, destfile = paste0(path_temp, "Art12_2018.gpkg"))

# Get Article 12 (Public!) - 2012
# Source: https://sdi.eea.europa.eu/catalogue/srv/eng/catalog.search#/metadata/7c2dd14f-60b6-4009-aca8-5d20300479a9
url <- "https://sdi.eea.europa.eu/webdav/datastore/public/eea_v_3035_10_mio_art12-2008-2012_p_2008-2012_v01_r01/Art12-2008-2012_GPKG/ART12_birds_distribution_reporting_0812_breeding_public.gpkg"
curl::curl_download(url, destfile = paste0(path_temp, "Art12_2012.gpkg"))

# Codes for Art 17
url <- "https://www.eea.europa.eu/data-and-maps/data/article-17-database-habitats-directive-92-43-eec-2/article-17-2020-dataset/article-17-2020-data-csv-format/at_download/file"
curl::curl_download(url, destfile = paste0(path_temp, "Article17_2020_dataset_csv.zip"))

# --- #
# EU Corine Accounting layer
## https://sdi.eea.europa.eu/catalogue/srv/eng/catalog.search#/metadata/a55d9224-a326-4cb1-9b9c-3a324520341a
# --> Downloaded manually to CLC folder
assertthat::assert_that(
  dir.exists(paste0(path_temp, "CLC_account"))
)
