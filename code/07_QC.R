# # # # # # # # # # # # # # # #
# Quality control checks for the outputs of 04_calculate_MTA.R
# Run from the project root after the MTA calculations.
#
# Author:
# Martin Jung (IIASA) - 2024-2026
# # # # # # # # # # # # # # # #

# Packages
library(cli)
library(assertthat)
library(dplyr)

path_output <- "export/"
version <- "02"

expected_files <- c(
  paste0("MTA_loglinear_EU_All_All_", version, ".csv"),
  paste0("MTA_loglinear_EU_Art12_species_", version, ".csv"),
  paste0("MTA_loglinear_EU_Art17_habitats_", version, ".csv"),
  paste0("MTA_loglinear_EU_Art17_species_", version, ".csv"),
  paste0("MTA_loglinear_MS_All_All_", version, ".csv"),
  paste0("MTA_loglinear_MS_Art12_species_", version, ".csv"),
  paste0("MTA_loglinear_MS_Art17_habitats_", version, ".csv"),
  paste0("MTA_loglinear_MS_Art17_species_", version, ".csv")
)
expected_paths <- file.path(path_output, expected_files)

qc_check <- function(condition, description) {
  if (!isTRUE(condition)) {
    cli::cli_abort("QC check failed: {description}")
  }
  cli::cli_alert_success("{description}")
}

cli::cli_alert_info("Checking MTA outputs in `{path_output}`.")
qc_check(all(file.exists(expected_paths)), "All expected MTA output files exist")

mta_tables <- lapply(expected_paths, read.csv, stringsAsFactors = FALSE)
names(mta_tables) <- expected_files
qc_check(all(vapply(mta_tables, nrow, integer(1)) > 0L),
         "All MTA output files contain rows")

required_columns <- c("option", "ctype", "category", "year", "y", "sd",
                      "variable", "scale", "art")
qc_check(all(vapply(mta_tables, function(x) all(required_columns %in% names(x)),
                    logical(1))),
         "All MTA output files contain the required columns")

mta <- dplyr::bind_rows(mta_tables)
expected_ctypes <- c("cdda", "combined", "natura2000")
expected_years <- c("<2000", "2000-2006", "2006-2012", "2012-2018", "2018-2024")
qc_check(all(vapply(mta_tables, function(x) setequal(unique(x$ctype), expected_ctypes),
                    logical(1))),
         "Every MTA file contains all expected conservation-area types")
qc_check(all(vapply(mta_tables, function(x) setequal(unique(x$year), expected_years),
                    logical(1))),
         "Every MTA file contains all expected reporting periods")

expected_metadata <- data.frame(
  scale = rep(c("EU", "MS"), each = 4),
  art = rep(c("Art_All", "Art_12", "Art_17", "Art_17"), times = 2),
  variable = rep(c("All", "species", "habitats", "species"), times = 2)
)
metadata_matches <- vapply(seq_along(mta_tables), function(i) {
  table <- mta_tables[[i]]
  all(unique(table$scale) == expected_metadata$scale[i],
      unique(table$art) == expected_metadata$art[i],
      unique(table$variable) == expected_metadata$variable[i])
}, logical(1))
qc_check(all(metadata_matches), "MTA metadata matches each output filename")

qc_check(!anyNA(mta), "MTA outputs contain no missing values")
qc_check(is.numeric(mta$y) && all(is.finite(mta$y)) &&
           all(mta$y >= 0 & mta$y <= 1),
         "All MTA values are within [0, 1]")
qc_check(is.numeric(mta$sd) && all(is.finite(mta$sd)) &&
           all(mta$sd >= 0 & mta$sd <= 1),
         "All MTA standard deviations are within [0, 1]")

group_columns <- c("option", "ctype", "category", "year", "variable", "scale", "art")
qc_check(anyDuplicated(mta[group_columns]) == 0L,
         "MTA output groups are unique")

feature_summary_file <- file.path(path_output, "Overall_statistics.csv")
qc_check(file.exists(feature_summary_file), "The feature-count summary exists")
feature_summary <- read.csv(feature_summary_file, stringsAsFactors = FALSE)
qc_check(all(c("category", "nr_features") %in% names(feature_summary)),
         "The feature-count summary contains the required columns")
species_counts <- feature_summary[which(feature_summary$category == "Species"), , drop = FALSE]
qc_check(nrow(species_counts) > 0L && !anyNA(species_counts$nr_features) &&
           all(species_counts$nr_features > 10),
         "Every species dataset contains more than 10 unique species")

comparison_columns <- c("option", "category", "year", "variable", "scale", "art")
combined <- mta[mta$ctype == "combined", c(comparison_columns, "y"), drop = FALSE]
names(combined)[names(combined) == "y"] <- "combined_y"

check_combined_ge_component <- function(ctype, label) {
  component <- mta[mta$ctype == ctype, c(comparison_columns, "y"), drop = FALSE]
  names(component)[names(component) == "y"] <- "component_y"
  comparison <- merge(combined, component, by = comparison_columns,
                      all.x = TRUE, sort = FALSE)

  qc_check(nrow(comparison) == nrow(combined) && !anyNA(comparison$component_y),
           paste("Combined MTA rows have matching", label, "rows"))
  qc_check(all(comparison$combined_y + 1e-10 >= comparison$component_y),
           paste("Combined MTA is not below", label, "MTA"))
}

check_combined_ge_component("natura2000", "N2k")
check_combined_ge_component("cdda", "CDDA")

cli::cli_alert_info("MTA QC passed for {length(expected_paths)} files and {nrow(mta)} rows.")
