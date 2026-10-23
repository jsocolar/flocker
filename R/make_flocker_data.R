# This file contains the function `make_flocker_data`, which calls one of the
# following (also contained in this file):
# make_flocker_data_static for a single-season model,
# make_flocker_data_dynamic for a multi-season model, or
# make_flocker_data_augmented for a data-augmented model, or
# make_flocker_data_twolevel_single for a two-level single-season model.


##### make_flocker_data ####
#' Format data for occupancy model with \code{flock()}.
#' @param obs If \code{type = "single"} or \code{type = "twolevel_single"},
#'  an I x J matrix-like object where closure is assumed across rows and columns
#'  are repeated sampling events.
#'    If \code{type = "multi"}, an I x J x K array where rows are sites or 
#'  species-sites, columns are repeated sampling events, and slices along the 
#'  third dimension are seasons. Allowable values are 1 (detection), 0 (no 
#'  detection), and NA (no sampling event).
#'     If \code{type = "augmented"}, an L x J x K array where rows L are sites, 
#'  columns J are repeat sampling events, and slices K are species. Every
#'  supplied species must have at least one detection; never-observed
#'  pseudospecies are added through \code{n_aug}. Species names may be supplied
#'  as names on the third dimension. If omitted, observed species are named
#'  `species1`, `species2`, and so on. Names matching
#'  `pseudospecies[digits]` are reserved for augmented pseudospecies.
#'     The data must be packed so that, for a given unit (site, site-species, 
#'  site-timestep, site-species-timestep) all realized visits come before any 
#'  missing visits (NAs are trailing within their rows).
#' @param unit_covs If \code{type = "single"} or
#'   \code{type = "twolevel_single"}, a dataframe of covariates for each
#'   closure unit that are constant across repeated sampling events within
#'   units. For \code{type = "twolevel_single"}, this must include the factor
#'   column named by \code{level2_group}.
#'   If \code{type = "multi"}, a list of such dataframes, one per timestep. All 
#' dataframes must have identical column names and types, and all
#' dataframes must have I rows.
#'   If \code{type = "augmented"}, a dataframe of covariates for each site that
#' are constant across repeated sampling events within sites (no dependence on
#' species is allowed).
#' @param event_covs If \code{type = "single"} or
#'   \code{type = "twolevel_single"}, a named list of I x J matrices, each one
#'   corresponding to a covariate that varies across repeated sampling events
#'   within closure-units.
#'   If \code{type = "multi"}, a named list of I x J x K arrays, each one 
#' corresponding to a covariate that varies across repeated sampling events 
#' within closure-units.
#'   If \code{type = "augmented"}, a named list of L x J matrices, each one
#' corresponding to a covariate that varies across repeated sampling events
#' within sites (no dependence on species is allowed).
#' @param type The type of occupancy model desired. Options are:
#'    \code{"single"} for a single-season model,
#'    \code{"multi"} for a multi-season (dynamic) model, or
#'    \code{"augmented"} for a single-season multi-species model with 
#'    data-augmentation for never-observed pseudospecies, or
#'    \code{"twolevel_single"} for a two-level single-season model.
#' @param n_aug Number of pseudo-species to augment. Only applicable if 
#'    \code{type = "augmented"}.
#' @param level2_group The name of a factor column in \code{unit_covs} and, if
#'    supplied, \code{level2_covs}, identifying the level-two group. Only
#'    applicable if \code{type = "twolevel_single"}.
#' @param level2_covs An optional dataframe with one row per level-two group
#'    and no unused levels in its grouping factor. Row order is arbitrary;
#'    factor-level order determines the internal group order, and every group
#'    must be represented in \code{unit_covs}. Only applicable if
#'    \code{type = "twolevel_single"}.
#' @param quiet Hide progress bars and informational messages?
#' @param newdata_checks If TRUE, turn off checks that must pass in order
#' to use the data for model fitting, but not in other contexts (e.g. making
#' predictions or assessing log-likelihoods over new data).
#' @return A flocker_data list that can be passed as data to \code{flock()}.
#' @details For one-level single-season data with event covariates and for
#' multiseason, two-level, and augmented data, the object
#' contains \code{data} and \code{metadata}. Observed visits are packed
#' contiguously within units, with unit covariates repeated at their actual
#' values on every visit. Consequently, brms centers unit-formula predictors
#' and constructs spline bases from covariates expanded to the visit level.
#' For multiseason data, units are additionally packed within series. Leading
#' and internal unvisited seasons receive one placeholder row, retaining their
#' unit covariates. Missing detection covariates on these rows are filled with
#' the mean over observed visits for numeric covariates, or the first factor
#' level for factors. Supplied values are preserved. Placeholder covariates
#' participate in brms preprocessing, but have no detection likelihood.
#' Trailing unvisited seasons are omitted and reconstructed as NA.
#' Two-level models pack units within groups using the same local block
#' structure. Group covariates are repeated across the group's visit rows and
#' therefore also participate in visit-level preprocessing. Units without
#' observations are not allowed in two-level models. Augmented data use the
#' same packing, with species as groups and site-species combinations as units.
#' Added pseudospecies share the observed species' sampling mask and covariates;
#' metadata reconstruct the site-by-visit-by-species arrays and species names.
#' Setting \code{newdata_checks = TRUE} marks this representation as
#' \code{new_data}, which cannot be used for fitting.
#' @export
#' @examples
#' sfd <- simulate_flocker_data()
#' make_flocker_data(
#'  sfd$obs, 
#'  sfd$unit_covs,
#'  sfd$event_covs
#' )
make_flocker_data <- function(obs, unit_covs = NULL, event_covs = NULL,
                              type = "single", n_aug = NULL,
                              quiet = FALSE, newdata_checks = FALSE,
                              level2_group = NULL, level2_covs = NULL) {
  standard_mfd_checks(obs, unit_covs, event_covs, type, n_aug, quiet, newdata_checks,
                      level2_group, level2_covs)

  if (!quiet) {
    if (type == "single") {
      message(paste0("Formatting data for a single-season occupancy model. For ",
                     "details, see make_flocker_data_static. All warnings and ",
                     "error messages should be interpreted in the context of ",
                     "make_flocker_data_static"))
    } else if (type == "multi") {
      message(paste0("Formatting data for a multi-season occupancy model. For ",
                     "details, see make_flocker_data_dynamic.  All warnings and ",
                     "error messages should be interpreted in the context of ",
                     "make_flocker_data_dynamic"))
    } else if (type == "augmented") {
      message(paste0("Formatting data for a single-season multispecies occupancy ", 
                     "model with data augmentation for never-observed species. For ",
                     "details, see make_flocker_data_augmented.  All warnings and ",
                     "error messages should be interpreted in the context of ",
                     "make_flocker_data_augmented"))
    } else if (type == "twolevel_single") {
      message(paste0("Formatting data for a two-level single-season occupancy ",
                     "model. For details, see make_flocker_data_twolevel_single. ",
                     "All warnings and error messages should be interpreted in ",
                     "the context of make_flocker_data_twolevel_single"))
    }
  }
  
  if (type == "single") {
    out <- make_flocker_data_static(
      obs, unit_covs, event_covs, quiet, newdata_checks
      )
  } else if (type == "multi") {
    out <- make_flocker_data_dynamic(
      obs, unit_covs, event_covs, quiet, newdata_checks
      )
  } else if (type == "augmented") {
    out <- make_flocker_data_augmented(
      obs, n_aug, unit_covs, event_covs, quiet, newdata_checks)
  } else if (type == "twolevel_single") {
    out <- make_flocker_data_twolevel_single(
      obs, unit_covs, event_covs, level2_group, level2_covs,
      quiet, newdata_checks)
  }
  
  out
}


##### make_flocker_data_static #####

#' Format data for single-season occupancy model, to be passed to 
#'  \code{flock()}.
#' @param obs An I x J matrix-like object where closure is assumed across rows 
#'  and columns are repeated sampling events. Allowable values are 1 (detection), 
#'  0 (no detection), and NA (no sampling event).
#   The data must be formatted so that all NAs are trailing within their rows.
#' @param unit_covs A dataframe of covariates for each unit that are constant 
#'            across repeated sampling events within closure-units.
#' @param event_covs A named list of I x J matrices, each one corresponding to a covariate
#' that varies across repeated sampling events within closure-units
#' @param quiet Hide progress bars and informational messages?
#' @param newdata_checks If TRUE, turn off checks that must pass in order
#' to use the data for model fitting, but not in other contexts (e.g. making
#' predictions or assessing log-likelihoods over new data).
#' @return A flocker_data list that can be passed as data to \code{flock()}.
#' @export
#' @examples
#' sfd <- simulate_flocker_data()
#' make_flocker_data_static(
#'  sfd$obs, 
#'  sfd$unit_covs,
#'  sfd$event_covs
#' )
make_flocker_data_static <- function(
    obs, unit_covs = NULL, event_covs = NULL, 
    quiet = FALSE, newdata_checks = FALSE
    ) {
  standard_mfd_checks(obs, unit_covs, event_covs, "single", NULL, quiet, newdata_checks)
  
  n_unit <- nrow(obs)
  n_rep <- ncol(obs)

  if (is.null(event_covs)) {
    n_trial <- rowSums(!is.na(obs))
    n_suc <- rowSums(obs, na.rm = T)
    flocker_data <- data.frame(ff_n_suc = n_suc, ff_n_trial = n_trial)
    if (!is.null(unit_covs)) {
      flocker_data <- cbind(flocker_data, unit_covs)
    }
    out <- list(data = flocker_data, metadata = list(
      type = "single_C", schema_version = 2L,
      data_role = if (newdata_checks) "new_data" else "model_data",
      flocker_version = flocker_version(), n_rep = n_rep,
      original_dim = dim(obs), original_dimnames = dimnames(obs),
      unit_covs = names(unit_covs), event_covs = names(event_covs)
    ))
  } else {
    # Pack actual events contiguously within each unit. Unit covariates retain
    # their real values on every visit, including when used for detection.
    obs <- as.matrix(obs)
    n_event <- as.integer(rowSums(!is.na(obs)))
    packed_to_original <- as.vector(t(matrix(seq_along(obs), nrow = n_unit)))
    packed_to_original <- packed_to_original[!is.na(obs[packed_to_original])]
    event_unit <- rep(seq_len(n_unit), n_event)
    flocker_data <- data.frame(ff_y = as.integer(obs[packed_to_original]))
    if (!is.null(unit_covs)) {
      flocker_data <- cbind(flocker_data, unit_covs[event_unit, , drop = FALSE])
    }
    for (name in names(event_covs)) {
      flocker_data[[name]] <- as.matrix(event_covs[[name]])[packed_to_original]
    }
    unit_start <- cumsum(c(1L, head(n_event, -1L)))
    n_data <- nrow(flocker_data)
    flocker_data$ff_n_unit <- c(n_unit, rep(-99L, n_data - 1L))
    flocker_data$ff_unit_start <- c(unit_start, rep(-99L, n_data - n_unit))
    flocker_data$ff_Q <- as.integer(rowSums(obs, na.rm = TRUE) > 0)[event_unit]
    rownames(flocker_data) <- NULL
    original_to_packed <- array(NA_integer_, dim(obs), dimnames(obs))
    original_to_packed[packed_to_original] <- seq_along(packed_to_original)
    out <- list(
      data = flocker_data,
      metadata = list(
        type = "single", schema_version = 2L,
        data_role = if (newdata_checks) "new_data" else "model_data",
        flocker_version = flocker_version(),
        original_dim = dim(obs), original_dimnames = dimnames(obs),
        packed_to_original = packed_to_original,
        original_to_packed = original_to_packed,
        event_unit = event_unit,
        unit_start = unit_start,
        n_event = n_event,
        unit_covs = names(unit_covs), event_covs = names(event_covs)
      )
    )
    class(out) <- c("list", "flocker_data")
    return(out)
  }
  class(out) <- c("list", "flocker_data")
  out
}


##### make_flocker_data_dynamic #####

#' Format data for dynamic (multi-season) occupancy model, to be passed to 
#'  \code{flock()}.
#' @param obs An I x J x K array where closure is assumed across rows, columns 
#'  are repeated sampling events, and slices along the third dimension are 
#'  seasons. Allowable values are 1 (detection), 0 (no detection), and NA (no 
#'  sampling event).
#'  The data must be formatted so that all NAs are trailing within their rows
#'  across repeat visits, but not necessarily across seasons. 
#' @param unit_covs  A list of dataframes (one per season) of covariates for 
#'    each closure-unit that are constant across repeated sampling events within 
#'    units. All dataframes must have identical column names and types, and all
#'    dataframes must have I rows.
#' @param event_covs A named list of I x J x K arrays, each one corresponding to 
#' a covariate that varies across repeated sampling events within closure-units
#' @param quiet Hide progress bars and informational messages?
#' @param newdata_checks If TRUE, turn off checks that must pass in order
#' to use the data for model fitting, but not in other contexts (e.g. making
#' predictions or assessing log-likelihoods over new data).
#' @return A flocker_data list that can be passed as data to \code{flock()}.
#' @export
make_flocker_data_dynamic <- function(obs, unit_covs = NULL, event_covs = NULL,
                                      quiet = FALSE, newdata_checks = FALSE) {
  
  standard_mfd_checks(obs, unit_covs, event_covs, "multi", NULL, quiet, newdata_checks)
  
  n_series <- nrow(obs)
  n_visit <- ncol(obs)
  n_season <- nslice(obs)
  first_visit <- matrix(obs[, 1, , drop = FALSE], n_series, n_season)
  n_season_observed <- apply(first_visit, 1, max_position_not_na)
  unit_grid <- matrix(seq_len(n_series * n_season), n_series, n_season)
  retained <- as.vector(col(unit_grid) <= n_season_observed)
  retained_units <- which(retained)
  # Matrix rows are conceptual series-season units; columns are visits.
  original_indices <- matrix(aperm(array(seq_along(obs), dim(obs)), c(1, 3, 2)),
                             n_series * n_season, n_visit)
  unit_obs <- matrix(obs[as.vector(original_indices)], n_series * n_season, n_visit)
  unit_obs <- unit_obs[retained, , drop = FALSE]
  covariates <- if (is.null(unit_covs)) NULL else do.call(rbind, unit_covs)[retained, , drop = FALSE]
  event_values <- lapply(event_covs, function(x) {
    values <- x[as.vector(original_indices[retained, , drop = FALSE])]
    dim(values) <- dim(unit_obs)
    values
  })
  packed <- pack_nested_blocks(unit_obs, covariates, event_values,
                                rep(seq_len(n_series), n_season)[retained], outer = "series")
  data <- packed$data

  packed_to_original <- original_indices[retained, , drop = FALSE][packed$packed_to_original]
  original_to_packed <- array(NA_integer_, dim(obs), dimnames(obs))
  original_to_packed[as.vector(original_indices[retained, , drop = FALSE])] <-
    as.vector(packed$original_to_packed)
  unit_to_packed <- matrix(NA_integer_, n_series, n_season,
                           dimnames = dimnames(obs)[c(1, 3)])
  unit_order <- retained_units[packed$unit_order]
  unit_to_packed[retained_units] <- packed$unit_to_packed
  out <- list(data = data, metadata = list(
    type = "multi", schema_version = 2L,
    data_role = if (newdata_checks) "new_data" else "model_data",
    flocker_version = flocker_version(),
    original_dim = dim(obs), original_dimnames = dimnames(obs),
    packed_to_original = packed_to_original,
    original_to_packed = original_to_packed, unit_to_packed = unit_to_packed,
    unit_order = unit_order, series_start = packed$group_start,
    unit_start = packed$unit_start, n_event = packed$n_event,
    event_unit = packed$row_unit, placeholder = packed$placeholder,
    unit_covs = names(unit_covs[[1]]), event_covs = names(event_covs)
  ))
  class(out) <- c("list", "flocker_data")
  out
}

##### make_flocker_data_augmented #####

#' Format data for data-augmented occupancy model, to be passed to 
#'  \code{flock()}.
#' @param obs An I x J x K array where rows I are sites, columns J are 
#'  repeat sampling events, and slices K are species. Allowable values are 1 
#'  (detection), 0 (no detection), and NA (no sampling event). Every supplied
#'  species must have at least one detection; never-observed pseudospecies are
#'  added through \code{n_aug}. Species names may be supplied as names on the
#'  third dimension. If omitted, observed species are named `species1`,
#'  `species2`, and so on. Names matching `pseudospecies[digits]` are
#'  reserved for augmented pseudospecies.
#'   The data must be formatted so that all NAs are trailing within their rows.
#' @param n_aug Number of pseudospecies to augment
#' @param site_covs A dataframe of covariates for each site that are constant 
#'            across repeated sampling events.
#' @param event_covs A named list of I x J matrices, each one corresponding to a 
#' covariate that varies across repeated sampling events within sites
#' @param quiet Hide progress bars and informational messages?
#' @param newdata_checks If TRUE, turn off checks that must pass in order
#' to use the data for model fitting, but not in other contexts (e.g. making
#' predictions or assessing log-likelihoods over new data).
#' @return A flocker_data list that can be passed as data to \code{flocker()}.
#' @export
make_flocker_data_augmented <- function(obs, n_aug, site_covs = NULL, 
                                        event_covs = NULL, quiet = FALSE, 
                                        newdata_checks = FALSE) {
  standard_mfd_checks(obs, site_covs, event_covs, "augmented", n_aug, quiet, newdata_checks)
  n_observed_species <- dim(obs)[3]
  observed_species_names <- dimnames(obs)[[3]]
  if (is.null(observed_species_names)) {
    observed_species_names <- paste0("species", seq_len(n_observed_species))
  }
  group_names <- c(
    observed_species_names,
    paste0("pseudospecies", seq_len(n_aug))
  )
  input_dimnames <- dimnames(obs)
  obs1 <- matrix(obs[, , 1, drop = FALSE], dim(obs)[1], dim(obs)[2])
  n_rep <- ncol(obs1)
  n_sp <- dim(obs)[3] + n_aug
  n_site <- dim(obs)[1]
  aug_slice <- obs1
  aug_slice[!is.na(aug_slice)] <- 0
  dim(aug_slice) <- c(dim(aug_slice), 1)
  for (i in 1:n_aug) {
    obs <- abind::abind(obs, aug_slice, along = 3)
  }
  obs <- expand_array_3D(obs)
  ff_species <- factor(rep(seq_len(n_sp), each = n_site), levels = seq_len(n_sp))
  ff_site <- rep(seq_len(n_site), n_sp)
  unit_covs <- data.frame(ff_species = ff_species)
  if (!is.null(site_covs)) {
    unit_covs <- cbind(unit_covs, stack_matrix(site_covs, n_sp))
  }
  level2_covs <- data.frame(
    ff_species = factor(seq_len(n_sp), levels = seq_len(n_sp))
  )
  event_covs2 <- NULL
  if (!is.null(event_covs)) {
    event_covs2 <- lapply(event_covs, function(x){stack_matrix(x, n_sp)})
  }
  
  out <- make_flocker_data_twolevel_single_(
    obs = obs,
    unit_covs = unit_covs,
    event_covs = event_covs2,
    level2_covs = level2_covs,
    level2_group = "ff_species", newdata_checks = newdata_checks
  )
  metadata <- out$metadata
  metadata$type <- "augmented"
  metadata$n_sp <- n_sp
  metadata$level2_group_names <- group_names
  metadata$unit_site <- ff_site[metadata$unit_order]
  metadata$original_dim <- as.integer(c(n_site, n_rep, n_sp))
  metadata$original_dimnames <- list(input_dimnames[[1]], input_dimnames[[2]], group_names)
  metadata$original_to_packed <- aperm(
    array(metadata$original_to_packed, c(n_site, n_sp, n_rep)), c(1, 3, 2))
  metadata$unit_to_packed <- matrix(metadata$unit_to_packed, n_site, n_sp)
  metadata$unit_covs <- names(site_covs)
  metadata$event_covs <- names(event_covs)
  out$metadata <- metadata
  
  class(out) <- c("list", "flocker_data")
  out
}

##### make_flocker_data_twolevel_single #####

#' Format data for two-level single-season occupancy model, to be passed to
#' \code{flock()}.
#' @inheritParams make_flocker_data
#' @param level2_group The name of a factor column in \code{unit_covs} and, if
#'   supplied, \code{level2_covs}, identifying the level-two group.
#' @param level2_covs An optional dataframe with one row per level-two group
#'   and no unused levels in its grouping factor. Row order is arbitrary;
#'   factor-level order determines the internal group order, and every group
#'   must be represented in \code{unit_covs}. If omitted, a minimal table is
#'   constructed from the groups represented in \code{unit_covs}.
#' @return A flocker_data list that can be passed as data to \code{flock()}.
#' @export
make_flocker_data_twolevel_single <- function(
    obs, unit_covs, event_covs = NULL, level2_group, level2_covs = NULL,
    quiet = FALSE, newdata_checks = FALSE
    ) {
  standard_mfd_checks(obs, unit_covs, event_covs, "twolevel_single", NULL,
                      quiet, newdata_checks, level2_group, level2_covs)

  make_flocker_data_twolevel_single_(
    obs, unit_covs, event_covs, level2_group, level2_covs, newdata_checks
  )
}

# Shared packing workhorse; each public formatter validates its own inputs.
make_flocker_data_twolevel_single_ <- function(
    obs, unit_covs, event_covs, level2_group, level2_covs,
    newdata_checks = FALSE
    ) {
  obs <- as.matrix(obs)
  event_covs <- lapply(event_covs, function(x) if (is.data.frame(x)) as.matrix(x) else x)
  if (is.null(level2_covs)) {
    group_levels <- levels(droplevels(unit_covs[[level2_group]]))
    level2_covs <- stats::setNames(data.frame(factor(group_levels, levels = group_levels)),
                                  level2_group)
  }
  group_levels <- levels(level2_covs[[level2_group]])
  level2_covs <- level2_covs[match(group_levels, as.character(level2_covs[[level2_group]])), , drop = FALSE]
  group_id <- match(as.character(unit_covs[[level2_group]]), group_levels)
  group_covs <- setdiff(names(level2_covs), level2_group)
  covariates <- cbind(unit_covs, level2_covs[group_id, group_covs, drop = FALSE])
  packed <- pack_nested_blocks(obs, covariates, event_covs, group_id)
  data <- packed$data
  n_group <- length(group_levels)
  group_Q <- vapply(seq_len(n_group), function(g) {
    as.integer(any(obs[group_id == g, , drop = FALSE] == 1, na.rm = TRUE))
  }, integer(1))
  data$ff_group_Q <- group_Q[group_id[packed$unit_order[packed$row_unit]]]
  out <- list(data = data, metadata = list(
    type = "twolevel_single", schema_version = 2L,
    data_role = if (newdata_checks) "new_data" else "model_data",
    flocker_version = flocker_version(), original_dim = dim(obs),
    original_dimnames = dimnames(obs), original_to_packed = packed$original_to_packed,
    unit_to_packed = packed$unit_to_packed, unit_order = packed$unit_order,
    group_start = packed$group_start, unit_start = packed$unit_start,
    unit_group = group_id, group_representatives = match(seq_len(n_group), group_id),
    level2_group = level2_group, level2_group_names = group_levels,
    level2_covs = group_covs, unit_covs = names(unit_covs), event_covs = names(event_covs)
  ))
  class(out) <- c("list", "flocker_data")
  out
}


##### mfd input checking #####
#' input checking for make_flocker_data
#' @inheritParams make_flocker_data
standard_mfd_checks <- function(
    obs, unit_covs, event_covs, type, n_aug, quiet, newdata_checks,
    level2_group = NULL, level2_covs = NULL
) {
  
  # unique() on a data.frame returns unique rows, not values; flatten first so
  # matrix-like data.frames pass the same value check as matrices and arrays.
  unique_y <- unique(as.vector(if (is.data.frame(obs)) as.matrix(obs) else obs))
  assertthat::assert_that(
    all(unique_y %in% c(0, 1, NA)),
    msg = "obs may only contain the values 0, 1, or NA"
  )
  
  assertthat::assert_that(
    is.null(event_covs) | is_named_list(event_covs),
    msg = "event_covs must be a named list or NULL"
    )
  
  assertthat::assert_that(
    type %in% flocker_data_input_types(),
    msg = paste0("Invalid type argument. Type given as '", type, "' but must ",
                 "be one of the following: ", 
                 paste(flocker_data_input_types(), collapse = ", "))
  )
  
  assertthat::assert_that(
    !any(names(unit_covs) %in% names(event_covs)) &
      !any(names(unit_covs[[1]]) %in% names(event_covs)),
    msg = "overlapping names detected between unit_covs and event_covs"
  )
  
  for (i in seq_along(flocker_reserved())) {
    assertthat::assert_that(
      !any(grepl(flocker_reserved()[i], names(unit_covs))),
      msg = paste0("names of unit_covs include a reserved string matching ",
                   "the following regular expression: ", 
                   flocker_reserved()[i])
    )
    assertthat::assert_that(
      !any(grepl(flocker_reserved()[i], names(unit_covs[[1]]))),
      msg = paste0("names of unit_covs include a reserved string matching ",
                   "the following regular expression: ", 
                   flocker_reserved()[i])
    )
    assertthat::assert_that(
      !any(grepl(flocker_reserved()[i], names(event_covs))),
      msg = paste0("names of event_covs include a reserved string matching ",
                   "the following regular expression: ", 
                   flocker_reserved()[i])
    )
  }
  
  assertthat::assert_that(
    is_one_logical(quiet),
    msg = "quiet must be TRUE or FALSE"
  )
  assertthat::assert_that(
    is_one_logical(newdata_checks),
    msg = "newdata_checks must be TRUE or FALSE"
  )
  
  #### type specific checks ####
  if (!quiet & !is.null(n_aug) & (type != "augmented")) {
    warning(paste0("n_aug is set but will be ignored for type = '", type, "'."))
  }
  
  #### single checks ####
  if(type == "single") {
    n_unit <- nrow(obs)
    n_rep <- ncol(obs)
    
    assertthat::assert_that(
      length(dim(obs)) == 2,
      msg = "in a single-season model, obs must have exactly two dimensions"
    )
    
    assertthat::assert_that(
      !any(is.na(obs[ , 1])), 
      msg = paste0("obs has NAs in its first column; this is not allowed in ", 
                   "single-season models"
      )
    )
    
    assertthat::assert_that(
      newdata_checks | (n_rep >= 2), 
      msg = paste0(
        "obs must contain at least two columns unless being used for newdata ",
        "(see newdata_checks argument)."
      )
    )
    
    if (n_rep > 2) {
      for (j in 2:(n_rep - 1)) {
        the_nas <- is.na(obs[ , j])
        if (any(the_nas)) {
          the_nas2 <- which(the_nas)
          assertthat::assert_that(
            all(is.na(obs[the_nas2, j+1])),
            msg = "Some rows of obs have non-trailing NAs"
          )
        }
      }
    }
    assertthat::assert_that(
      !all(is.na(obs[ , n_rep])),
      msg = "The final column of obs contains only NAs."
    )
    if (!is.null(unit_covs)) {
      assertthat::assert_that(
        nrow(unit_covs) == nrow(obs),
        msg = "Different numbers of rows found for obs and unit_covs."
      )
      assertthat::assert_that(
        !any(is.na(unit_covs)),
        msg = "A unit covariate contains missing values."
      )
    }
    if (!is.null(event_covs)) {
      assertthat::assert_that(
        is_named_list(event_covs), 
        msg = "event_covs must be NULL or a named list with no duplicate names."
      )
      n_event_covs <- length(event_covs)
      missing_covs <- vector()
      for (ec in 1:n_event_covs) {
        assertthat::assert_that(
          all.equal(dim(event_covs[[ec]]), dim(obs)),
          msg = paste0(
            "Dimension mismatch found between obs and event_covs[[", ec, "]]."
          )
        )
        missing_covs <- unique(c(missing_covs, which(is.na(event_covs[[ec]]))))
      }
      if (length(missing_covs) > 0) {
        assertthat::assert_that(
          all(missing_covs %in% which(is.na(obs))),
          msg = paste0("An event covariate contains missing values ",
                       "at a position where the response is not missing.")
        )
      }
    }
  }
  
  #### multi checks ####
  if(type == "multi") {
    n_year <- nslice(obs) # nslice checks that obs is a 3-D array
    n_series <- nrow(obs)
    n_rep <- ncol(obs)
    n_total <- n_year*n_series*n_rep
    
    assertthat::assert_that(
      newdata_checks | (n_year > 1), 
      msg = paste0(
        "obs must contain at least two slices (seasons/years) or you must ",
        "format data explicitly as newdata (see newdata_checks argument)."
      )
    )
    assertthat::assert_that(
      newdata_checks | (n_rep > 1), 
      msg = paste0(
        "obs must contain at least two columns (repeat visits to at ",
        "least one unit), or you mus format data explicitly as newdata (see ",
        "newdata_checks argument)."
      )
    )
    
    # Check that no NAs are non-trailing across columns (i.e. reps within series-
    # years)
    for (k in seq(n_year)) {
      for (j in seq_len(n_rep - 1L)) {
        the_nas <- is.na(obs[ , j, k])
        if (any(the_nas)) {
          the_nas2 <- which(the_nas)
          assertthat::assert_that(
            all(is.na(obs[the_nas2, j+1, k])),
            msg = paste0("Some rows/slices of obs have non-trailing NAs ", 
                         "across columns."
            )
          )
        }
      }
    }
    
    # Check that the first and final reps contain at least one non-NA
    assertthat::assert_that(
      !all(is.na(obs[ , 1, ])),
      msg = "The first column (replicate visit) of obs contains only NAs."
    )
    assertthat::assert_that(
      (!all(is.na(obs[ , n_rep, ]))) | newdata_checks,
      msg = "The final column (replicate visit) of obs contains only NAs."
    )
    assertthat::assert_that(
      is.null(unit_covs) | is.list(unit_covs),
      msg = "unit_covs must be a list or NULL."
    )
    assertthat::assert_that(
      is.null(event_covs) | is_named_list(event_covs),
      msg = "event_covs must be a named list or NULL."
    )
    if (all(is.na(obs[ , , 1]))) {
      warning("The first slice (season) of obs contains only NAs")
    }  
    if (all(is.na(obs[ , , n_year]))) {
      warning("The final slice (season) of obs contains only NAs")
    }
    if (!is.null(unit_covs)) {
      assertthat::assert_that(
        length(unit_covs) == n_year,
        msg = "If provided, unit_covs must have length equal to dim(obs)[3]"
      )
      for (k in 1:n_year) {
        assertthat::assert_that(
          is.data.frame(unit_covs[[k]]),
          msg = "All elements of unit_covs must be dataframes."
        )
        assertthat::assert_that(
          identical(dim(unit_covs[[k]]), dim(unit_covs[[1]])),
          msg = "All elements of unit_covs must have identical dimensions."
        )
        assertthat::assert_that(
          identical(names(unit_covs[[k]]), names(unit_covs[[1]])),
          msg = "All elements of unit_covs must have identical column names."
        )
        assertthat::assert_that(
          !any(is.na(unit_covs[[k]])),
          msg = paste0("NA unit covariates are not allowed in dynamic models. ",
                       "It is safe to impute dummy values in the following ",
                       "circumstances.",
                       "Note, however, that imputing values can interfere ",
                       "with brms's default behavior of centering the ", 
                       "columns of the design matrix. To avoid nonintuitive ",
                       "prior specifications for the intercepts, impute the ", 
                       "mean value rather than any other choice of dummy.", 
                       "1) the model uses explicit initial occupancy ",
                       "probabilities, and a unit covariate is used only for ",
                       "initial occupancy and not for detection, ",
                       "colonization or extinction; impute values for years ",
                       "after the first. ",
                       "2) the model uses explicit initial occupancy ",
                       "probabilities, and a unit covariate is used only for ",
                       "colonization/extinction and not for initial ",
                       "occupancy or detection; impute values for the first ",
                       "year. ",
                       "3) a unit covariate is used only for detection; ",
                       "impute values for the first visit at units with no ",
                       "visits. ",
                       "4) a unit covariate for colonization or extinction ",
                       "is unavailable at a timestep with no observed data ",
                       "at the end of the timeseries, or a timestep that is ",
                       "part of a block of timesteps with no observed data ",
                       "reaching uninterrupted to the and of the timeseries, ",
                       "and inference on likely occupancy probabilties is ", 
                       "not desired at any of those timesteps; impute values ",
                       "for the trailing block of timesteps with no observations.")
        )
      }
      assertthat::assert_that(
        nrow(unit_covs[[1]]) == n_series,
        msg = "each element of unit_covs must have the same number of rows as obs"
      )
    }
    if (!is.null(event_covs)) {
      n_event_covs <- length(event_covs)
      missing_covs <- vector()
      for (ec in 1:n_event_covs) {
        assertthat::assert_that(
          all.equal(dim(event_covs[[ec]]), dim(obs)),
          msg = paste0("Dimension mismatch found between obs and event_covs[[", ec, "]].")
        )
        missing_covs <- unique(c(missing_covs, which(is.na(event_covs[[ec]]))))
      }
      assertthat::assert_that(
        all(missing_covs %in% which(is.na(obs))),
        msg = paste0("An event covariate contains missing values ",
                     "at a position where the response is not missing.")
      )
    }
    assertthat::assert_that(
      !is.null(event_covs),
      msg = paste0("Construction alert! The model contains no event covariates. ",
                   "This is fine, but for now please add a dummy event covariate.",
                   "You do not need to use this covariate in your model formula.")
    )
    
    
    n_year_obs <- apply(matrix(obs[, 1, , drop = FALSE], n_series, n_year),
                        1, max_position_not_na)
    assertthat::assert_that(
      !any(n_year_obs == 0),
      msg = paste0("at least one series (i.e. row; generally a site or a ",
                   "species-site) has no observations at any timestep")
    )
  }
  
  #### twolevel_single checks ####
  if(type == "twolevel_single"){
    assertthat::assert_that(
      length(dim(obs)) == 2,
      msg = "in a two-level single-season model, obs must have exactly two dimensions"
    )
    assertthat::assert_that(
      !is.null(unit_covs),
      msg = "unit_covs must be supplied for two-level single-season models"
    )
    assertthat::assert_that(
      is.data.frame(unit_covs),
      msg = "unit_covs must be a dataframe"
    )
    assertthat::assert_that(
      is.null(level2_covs) | is.data.frame(level2_covs),
      msg = "level2_covs must be NULL or a dataframe"
    )
    assertthat::assert_that(
      is.character(level2_group) & length(level2_group) == 1,
      msg = "level2_group must be a single column name"
    )
    assertthat::assert_that(
      level2_group %in% names(unit_covs),
      msg = "level2_group must name a column in unit_covs"
    )
    assertthat::assert_that(
      is.factor(unit_covs[[level2_group]]),
      msg = "level2_group must identify a factor column in unit_covs"
    )
    assertthat::assert_that(
      newdata_checks | anyDuplicated(unit_covs[[level2_group]]) > 0,
      msg = "At least one level-two group must contain more than one unit."
    )
    if (!is.null(level2_covs)) {
      assertthat::assert_that(
        level2_group %in% names(level2_covs),
        msg = "level2_group must name a column in level2_covs"
      )
      assertthat::assert_that(
        is.factor(level2_covs[[level2_group]]),
        msg = paste0(
          "level2_group must identify a factor column in level2_covs"
        )
      )
      duplicate_level2_covs <- setdiff(
        intersect(names(unit_covs), names(level2_covs)),
        level2_group
      )
      assertthat::assert_that(
        length(duplicate_level2_covs) == 0,
        msg = paste0(
          "level2_covs and unit_covs may only share the level2_group column. ",
          "Duplicate column name(s): ",
          paste(duplicate_level2_covs, collapse = ", ")
        )
      )
      assertthat::assert_that(
        !any(names(level2_covs) %in% names(event_covs)),
        msg = "overlapping names detected between level2_covs and event_covs"
      )
      for (reserved in flocker_reserved()) {
        assertthat::assert_that(
          !any(grepl(reserved, names(level2_covs))),
          msg = paste0(
            "names of level2_covs include a reserved string matching ",
            "the following regular expression: ", reserved
          )
        )
      }
      assertthat::assert_that(
        !anyDuplicated(level2_covs[[level2_group]]),
        msg = "level2_covs must contain no duplicate level2_group values"
      )
      assertthat::assert_that(
        setequal(
          as.character(level2_covs[[level2_group]]),
          levels(level2_covs[[level2_group]])
        ),
        msg = paste0(
          "level2_covs must contain exactly one row for each level2_group ",
          "factor level"
        )
      )
      assertthat::assert_that(
        all(
          as.character(unit_covs[[level2_group]]) %in%
            as.character(level2_covs[[level2_group]])
        ),
        msg = "all level2_group values in unit_covs must appear in level2_covs"
      )
      empty_groups <- setdiff(
        levels(level2_covs[[level2_group]]),
        as.character(unit_covs[[level2_group]])
      )
      if (length(empty_groups) > 0) {
        stop(
          paste0(
            "level2_covs contains groups with no corresponding units: ",
            paste(empty_groups, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      assertthat::assert_that(
        !any(is.na(level2_covs)),
        msg = "A level-two covariate contains missing values."
      )
    }
    assertthat::assert_that(
      nrow(unit_covs) == nrow(obs),
      msg = "Different numbers of rows found for obs and unit_covs."
    )
    assertthat::assert_that(
      !any(is.na(unit_covs)),
      msg = "A unit covariate contains missing values."
    )
    assertthat::assert_that(
      !any(is.na(obs[ , 1])), 
      msg = paste0("obs has NAs in its first column; this is not allowed in ", 
                   "two-level single-season models")
    )
    assertthat::assert_that(
      newdata_checks | (ncol(obs) >= 2), 
      msg = paste0(
        "obs must contain at least two columns unless being used for newdata ",
        "(see newdata_checks argument)."
      )
    )
    if (ncol(obs) > 2) {
      for (j in 2:(ncol(obs) - 1)) {
        the_nas <- is.na(obs[ , j])
        if (any(the_nas)) {
          the_nas2 <- which(the_nas)
          assertthat::assert_that(
            all(is.na(obs[the_nas2, j+1])),
            msg = "Some rows of obs have non-trailing NAs"
          )
        }
      }
    }
    assertthat::assert_that(
      !all(is.na(obs[ , ncol(obs)])),
      msg = "The final column of obs contains only NAs."
    )
    if (!is.null(event_covs)) {
      assertthat::assert_that(
        is_named_list(event_covs), 
        msg = "event_covs must be NULL or a named list with no duplicate names."
      )
      missing_covs <- vector()
      for (ec in seq_along(event_covs)) {
        assertthat::assert_that(
          all.equal(dim(event_covs[[ec]]), dim(obs)),
          msg = paste0(
            "Dimension mismatch found between obs and event_covs[[", ec, "]]."
          )
        )
        missing_covs <- unique(c(missing_covs, which(is.na(event_covs[[ec]]))))
      }
      if (length(missing_covs) > 0) {
        assertthat::assert_that(
          all(missing_covs %in% which(is.na(obs))),
          msg = paste0("An event covariate contains missing values ",
                       "at a position where the response is not missing.")
        )
      }
    }
  }
  
  #### augmented checks ####
  if(type == "augmented"){
    site_covs <- unit_covs
    obs1 <- obs[,,1]
    n_rep <- ncol(obs1)
    na_obs <- which(is.na(obs1))
    
    assertthat::assert_that(
      length(dim(obs)) == 3,
      msg = "obs must have exactly three dimensions."
    )
    assertthat::assert_that(
      dim(obs)[1] > 1,
      msg = "Augmented models require more than one site."
    )
    assertthat::assert_that(
      is_one_pos_int(n_aug),
      msg = "n_aug must be a positive integer"
    )
    species_names <- dimnames(obs)[[3]]
    assertthat::assert_that(
      is.null(species_names) |
        !any(grepl("^pseudospecies[0-9]+$", species_names)),
      msg = paste0(
        "Observed species names matching 'pseudospecies' followed by digits ",
        "are reserved for augmented pseudospecies"
      )
    )
    detected_species <- apply(obs == 1, 3, any, na.rm = TRUE)
    if (any(!detected_species)) {
      stop(
        paste0(
          "The original obs array for an augmented model must contain only ",
          "species with at least one detection. Add never-observed ",
          "pseudospecies through n_aug."
        ),
        call. = FALSE
      )
    }
    
    for (i in 2:dim(obs)[3]) {
      na_obs_i <- which(is.na(obs[,,i]))
      assertthat::assert_that(
        identical(na_obs, na_obs_i),
        msg = "Different species have different sampling events NA"
      )
    }
    
    assertthat::assert_that(
      n_rep >= 2,
      msg = "obs must contain at least two columns."
    )

    assertthat::assert_that(
      all(!is.na(obs1[ , 1])),
      msg = "Some sites have NAs on the first sampling event."
    )
    if (n_rep > 2) {
      for (j in 2:(n_rep - 1)) {
        the_nas <- is.na(obs1[ , j])
        if (any(the_nas)) {
          the_nas2 <- which(the_nas)
          assertthat::assert_that(all(is.na(obs1[the_nas2, j+1])),
                                  msg = "Some sites have non-trailing NA visits."
          )
        }
      }
    }
    assertthat::assert_that(!all(is.na(obs1[ , n_rep])),
                            msg = "The final repeat event contains only NAs."
    )
    
    if (!is.null(site_covs)) {
      assertthat::assert_that(
        nrow(site_covs) == nrow(obs1),
        msg = "Different numbers of rows found for obs and site_covs."
      )
      assertthat::assert_that(
        all(!is.na(site_covs)),
        msg = "A site covariate contains missing values."
      )
    }
    if (!is.null(event_covs)) {
      assertthat::assert_that(
        is_named_list(event_covs),
        msg = "event_covs must be NULL or a named list with unique names"
      )
      n_event_covs <- length(event_covs)
      missing_covs <- vector()
      for (ec in 1:n_event_covs) {
        assertthat::assert_that(
          identical(dim(event_covs[[ec]]), dim(obs1)),
          msg = paste0("Dimension mismatch found between obs and event_covs[[", ec, "]].")
        )
        missing_covs <- unique(c(missing_covs, which(is.na(event_covs[[ec]]))))
      }
      assertthat::assert_that(
        all(missing_covs %in% which(is.na(obs1))),
        msg = paste0("An event covariate contains missing values ",
                     "at a position where the response is not missing.")
      )
    }
  }
}
