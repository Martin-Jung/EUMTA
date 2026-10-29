# This script creates broad-group CORINE masks at 100 m and their fractional
# coverage on the European 10 km reference grid.

library(terra)
library(assertthat)
library(stringr)

# Check for local GDAL installation and gdal_calc.py availability
is_windows <- .Platform$OS.type == "windows"
if (is_windows) {
  osgeo4w <- Sys.getenv("OSGEO4W_ROOT", unset = "C:/OSGeo4W")
  gdal_python <- file.path(osgeo4w, "apps", "Python312", "python.exe")
  gdal_calc <- file.path(osgeo4w, "apps", "Python312", "Scripts", "gdal_calc.py")
  if (!file.exists(gdal_python) || !file.exists(gdal_calc)) {
    stop("OSGeo4W (QGIS) must be installed on Windows to run gdal_calc.py. ",
         "Expected files under ", osgeo4w, ".")
  }
} else {
  gdal_calc <- Sys.which("gdal_calc.py")
  if (!nzchar(gdal_calc)) {
    stop("gdal_calc.py must be installed and available on PATH.")
  }
}

path_temp <- "temp/"
dir.create(path_temp, showWarnings = FALSE, recursive = TRUE)
path_corine <- file.path(path_temp, "temporary_data/CLC_account")
assertthat::assert_that(dir.exists(path_corine))

# CORINE accounting filenames contain their timestep in the CLC####ACC prefix.
ll_corine <- list.files(path_corine, pattern = "\\.tif(f)?$",
                        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
ll_corine <- ll_corine[grepl("^CLC[0-9]{4}ACC.*\\.tif$", basename(ll_corine),
                             ignore.case = TRUE)]
corine_year <- stringr::str_match(basename(ll_corine), "^CLC([0-9]{4})ACC")[, 2]
df_corine <- data.frame(file = ll_corine, year = corine_year,
                        stringsAsFactors = FALSE)
assertthat::assert_that(nrow(df_corine) > 0, !anyNA(df_corine$year),
                        msg = "No CORINE accounting GeoTIFFs found in temporary_data/CLC_account.")
df_corine <- df_corine[order(as.integer(df_corine$year)), , drop = FALSE]

# GRID_CODE stores raster values; CLC_CODE identifies their broad group.
corine_legend <- read.csv2(file.path(path_corine, "CLC_legend.csv"),
                           stringsAsFactors = FALSE)
assertthat::assert_that(all(c("GRID_CODE", "CLC_CODE") %in% names(corine_legend)),
                        msg = "CLC_legend.csv must contain GRID_CODE and CLC_CODE columns.")
corine_legend$CLC_group <- substr(as.character(corine_legend$CLC_CODE), 1, 1)
group_labels <- c("artificial_surfaces", "agricultural_areas",
                  "forest", "wetlands", "water_bodies")
names(group_labels) <- as.character(1:5)

background_10km <- terra::rast(
  "../11.1.2_ReferenceGrids/ReferenceGrid_Europe_bin_10000m.tif"
)

# Process every timestep and select all raster codes belonging to each group.
for (i in seq_len(nrow(df_corine))) {
  yr <- df_corine$year[i]
  input_file <- df_corine$file[i]
  message("Processing CORINE year: ", yr)

  for (group in names(group_labels)) {
    class_codes <- unique(as.integer(
      corine_legend$GRID_CODE[corine_legend$CLC_group == group]
    ))
    assertthat::assert_that(length(class_codes) > 0,
                            msg = paste("No CORINE raster codes found for group", group))

    group_label <- group_labels[[group]]
    output_100m <- file.path(path_temp, paste0("CLC_", group_label, "_", yr, "_100m.tif"))
    output_10km <- file.path(path_temp, paste0("CLC_", group_label, "_", yr, "_10km.tif"))
    calc_expression <- paste0("(", paste0("A==", class_codes, collapse = " | "), ")")
    args <- c("-A", input_file, paste0("--outfile=", output_100m),
              paste0("--calc=", calc_expression), "--NoDataValue=255", "--type=Byte",
              "--hideNoData", "--overwrite", "--co=COMPRESS=DEFLATE", "--co=TILED=YES")

    if (is_windows) {
      status <- system2(gdal_python, args = shQuote(c(gdal_calc, args), type = "cmd"))
    } else {
      status <- system2(gdal_calc, args = shQuote(args))
    }
    if (status != 0L) {
      stop("gdal_calc.py failed for year ", yr, ", group ", group_label,
           " (exit status ", status, ")")
    }

    group_100m <- terra::rast(output_100m)
    group_10km <- terra::project(group_100m, background_10km, method = "average")
    group_10km <- terra::mask(group_10km, background_10km)
    terra::writeRaster(group_10km, output_10km, datatype = "FLT4S",
                       overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
    rm(group_100m, group_10km)
  }
}