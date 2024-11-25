
#' Target calculations for a given feature
#'
#' @description
#' This function calculates the targets
#'
#' @author Martin Jung
#' @keywords targets, indicators
calc_targets <- function(data = NULL,
                         current_range, potential_range,
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

  # Drop units if existing
  if(inherits(current_range, "units")) current_range <- units::drop_units(current_range)
  if(inherits(potential_range, "units")) potential_range <- units::drop_units(potential_range)

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
                      target_relative = (potential_range * default_target)/potential_range
                      )
  } else if(option == "loglinear"){
    # Log-linear targets
    tr <- prioritizr:::loglinear_interpolation(
      potential_range,
      1000, # lower_bound_amount,
      1,  # lower_bound_target
      250000, # upper_bound_amount,
      0.1 # upper_bound_target
    ) * potential_range

    # Cap to 1 million km2
    if(any(tr>1e6)) tr[tr>1e6] <- 1e6
    # Format output
    out <- data.frame(option = option,
                      target_absolute = tr,
                      target_relative = (tr / potential_range) )
  } else if(option == "extinctrisk"){
    # Resolution is in km!
    assertthat::assert_that(
      is.data.frame(data),utils::hasName(data, "code")
    )
    out <- data |>
      units::drop_units() |>
      dplyr::group_by(ctype, code) |>  # Per species
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

  ## Set correct units
  # Set to km2
  out$target_absolute <- units::set_units(out$target_absolute, "km2")

  # Return targets
  assertthat::assert_that(is.data.frame(out), !anyNA(out),
                          utils::hasName(out, "target_relative"),
                          utils::hasName(out, "target_absolute")
                          )
  return(out)
}

#' MTA function
#' @description
#' Small helper function to calculate the MTA for a given set of data.
#' @param target_relative A [`numeric`] value between 0 and 1 indicating the target
#' for a given feature.
#' @param proportion A [`numeric`] value between 0 and 1 indicating the proportion
#' of conserved area.
#' @param target_absolute A [`numeric`] alternative value describing an absolute target (Default: \code{NULL}).
#' @param conserved_absolute A [`numeric`] alternative value with the conserved amount (Default: \code{NULL}).
#' @returns A [`numeric`] value.
#' @keywords indicator
mta <- function(target_relative = NULL, proportion = NULL,
                target_absolute = NULL, conserved_absolute = NULL){
  assertthat::assert_that(
    is.null(target_relative) || is.numeric(target_relative),
    # For proportion
    is.null(proportion) || is.numeric(proportion)
  )

  if(!is.null(target_absolute)){
    if(inherits(target_absolute, "unit")) target_absolute <- units::drop_units(target_absolute)
    if(inherits(conserved_absolute, "unit")) conserved_absolute <- units::drop_units(conserved_absolute)

    return(
      mean(
        pmin(1, conserved_absolute / target_absolute)
      )
    )

  } else {
    assertthat::assert_that(
      all( dplyr::between(target_relative, 0, 1)),
      all( dplyr::between(proportion, 0, 1))
    )
    return(
      mean(
        1- (pmin(1, pmax(0, ( target_relative - proportion ) )) / target_relative )
      )
    )
  }
}
