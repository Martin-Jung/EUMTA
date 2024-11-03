
#' Target calculations for a given feature
#'
#' @description
#' This function calculates the targets
#'
#' @author Martin Jung
#' @keywords targets, indicators
calc_targets <- function(current_range, potential_range,
                         data = NULL,
                         option = "flat", default_target = 0.3){
  assertthat::assert_that(
    is.character(option),
    # Check numeric range
    is.null(current_range) || (is.numeric(current_range) && length(current_range)>=1),
    is.numeric(potential_range), length(potential_range)>=1
  )
  # Match target calculation option
  option <- match.arg(option, c("none", "flat", "loglinear", "extinctrisk"),
                      several.ok = FALSE)

  # Flat targets
  if(option == "flat"){
    assertthat::assert_that(
      is.numeric(default_target),
      dplyr::between(default_target,0,1)
    )
    message("Computing flat targets with: ", default_target)

    # Simply multiply the potential ranges with the flat target
    out <- data.frame(option = option,
                      target_absolute = potential_range * default_target,
                      target_relative = (potential_range * default_target)/potential_range)
  } else if(option == "loglinear"){
    # Log-linear targets
    tr <- prioritizr:::loglinear_interpolation(
      potential_range |> units::drop_units(),
      1000, # lower_bound_amount,
      0.9,  # lower_bound_target
      250000, # upper_bound_amount,
      0.2 # upper_bound_target
    ) * potential_range |>
      units::drop_units()

    # Cap to 1 million km2
    if(any(tr>1e6)) tr[tr>1e6] <- 1e6
    # Set to km2
    tr <- units::set_units(tr, "km2")

    # Format output
    out <- data.frame(option = option,
                      target_absolute = tr,
                      target_relative = (tr / potential_range) |> units::drop_units())

  } else if(option == "extinctrisk"){
    # Resolution is in km!
    assertthat::assert_that(
      is.data.frame(data),utils::hasName(data, "code"),
      units::deparse_unit(current_range) == "km2",
      anyDuplicated(data$code) == 0
    )
    out <- data |>
      units::drop_units() |>
      dplyr::group_by(code) |>  # Per species
      dplyr::reframe(
        size = sum(totalarea_km2),
        target_absolute = min( c( max( c(22000), 0.8 * size  ) ), 10^6 )
      ) |>
      ungroup() |>
      # Add relative target as well
      dplyr::mutate(option = option,
                    target_relative = ifelse(  (target_absolute / size)>=1, 1, (target_absolute / size) ) ) |>
      ungroup()
  }
  # Return targets
  assertthat::assert_that(is.data.frame(out), !anyNA(out),
                          utils::hasName(out, "target_relative"),
                          utils::hasName(out, "target_absolute") )
  return(out)
}
