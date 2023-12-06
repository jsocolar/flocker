# This file contains the function `make_flocker_data`, which calls one of the
# following (also contained in this file):
# make_flocker_data_static for a single-season model,
# make_flocker_data_dynamic for a multi-season model, or
# make_flocker_data_augmented for a data-augmented model.


##### make_flocker_data ####
#' Format data for occupancy model with \code{flock()}.
#' @param obs If \code{type = "single"}, an I x J matrix-like object where 
#'  closure is assumed across rows and columns are repeated sampling events. 
#'    If \code{type = "multi"}, an I x J x K array where rows are sites or 
#'  species-sites, columns are repeated sampling events, and slices along the 
#'  third dimension are seasons. Allowable values are 1 (detection), 0 (no 
#'  detection), and NA (no sampling event).
#'     If \code{type = "augmented"}, an L x J x K array where rows L are sites, 
#'  columns J are repeat sampling events, and slices K are species. 
#'     The data must be packed so that, for a given unit (site, site-species, 
#'  site-timestep, site-species-timestep) all realized visits come before any 
#'  missing visits (NAs are trailing within their rows).
#' @param unit_covs If \code{type = "single"} a dataframe of covariates for each 
#' closure-unit that are constant across repeated sampling events within units.
#'   If \code{type = "multi"}, a list of such dataframes, one per timestep. All 
#' dataframes must have identical column names and types, and all
#' dataframes must have I rows.
#'   If \code{type = "augmented"}, a dataframe of covariates for each site that
#' are constant across repeated sampling events within sites (no dependence on
#' species is allowed).
#' @param event_covs If \code{type = "single"}, a named list of I x J matrices, 
#' each one corresponding to a covariate that varies across repeated sampling 
#' events within closure-units.
#'   If \code{type = "multi"}, a named list of I x J x K arrays, each one 
#' corresponding to a covariate that varies across repeated sampling events 
#' within closure-units.
#'   If \code{type = "augmented"}, a named list of L x J matrices, each one
#' corresponding to a covariate that varies across repeated sampling events
#' within sites (no dependence on species is allowed).
#' @param type The type of occupancy model desired. Options are:
#'    \code{"single"} for a single_season model,
#'    \code{"multi"} for a multi-season (dynamic) model, or
#'    \code{"augmented"} for a single-season multi-species model with 
#'    data-augmentation for never-observed pseudospecies.
#' @param n_aug Number of pseudo-species to augment. Only applicable if 
#'    \code{type = "augmented"}.
#' @param quiet Hide progress bars and informational messages?
#' @return A flocker_data list that can be passed as data to \code{flock()}.
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
                              quiet = FALSE) {
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
    }
    if (!is.null(n_aug) & (type != "augmented")) {
      warning(paste0("n_aug is set but will be ignored for type = '", type, "'."))
    }
  }
  
  if (type == "single") {
    out <- make_flocker_data_static(obs, unit_covs, event_covs, quiet)
    out$unit_covs <- names(unit_covs)
    out$event_covs <- names(event_covs)
  } else if (type == "multi") {
    out <- make_flocker_data_dynamic(obs, unit_covs, event_covs, quiet)
    out$unit_covs <- names(unit_covs[[1]])
    out$event_covs <- names(event_covs)
  } else if (type == "augmented") {
    out <- make_flocker_data_augmented(obs, n_aug, unit_covs, event_covs, quiet)
    out$unit_covs <- names(unit_covs)
    out$event_covs <- names(event_covs)
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
#' @return A flocker_data list that can be passed as data to \code{flock()}.
#' @export
#' @examples
#' sfd <- simulate_flocker_data()
#' make_flocker_data_static(
#'  sfd$obs, 
#'  sfd$unit_covs,
#'  sfd$event_covs
#' )
make_flocker_data_static <- function(obs, unit_covs = NULL, event_covs = NULL, quiet = FALSE) {
  assertthat::assert_that(
    length(dim(obs)) == 2,
    msg = "in a single-season model, obs must have exactly two dimensions"
  )
  n_unit <- nrow(obs)
  n_rep <- ncol(obs)
  assertthat::assert_that(
    n_rep >= 2, 
    msg = "obs must contain at least two columns."
  )
  unique_y <- unique(obs)[!is.na(unique(obs))]
  assertthat::assert_that(
    all(unique_y %in% c(0, 1)),
    msg = "obs may only contain the values 0, 1, or NA"
  )
  assertthat::assert_that(
    !any(is.na(obs[ , 1])), 
    msg = paste0("obs has NAs in its first column; this is not allowed in ", 
                 "single-season models"
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
  if (is.null(event_covs)) {
    n_trial <- rowSums(!is.na(obs))
    n_suc <- rowSums(obs, na.rm = T)
    flocker_data <- data.frame(ff_n_suc = n_suc, ff_n_trial = n_trial)
    if (!is.null(unit_covs)) {
      flocker_data <- cbind(flocker_data, unit_covs)
    }
    out <- list(data = flocker_data, n_rep = n_rep, 
                type = "single_C")
  } else {
    flocker_data <- data.frame(ff_y = expand_matrix(obs))
    if (!is.null(unit_covs)) {
      unit_covs_stacked <- 
        do.call(rbind, replicate(n_rep, unit_covs, simplify=FALSE))
      flocker_data <- cbind(flocker_data, unit_covs_stacked)
    }
    event_covs <- as.data.frame(lapply(event_covs, expand_matrix))
    flocker_data <- cbind(flocker_data, event_covs)
    flocker_data$ff_n_unit <- c(nrow(obs), 
                             rep(-99, nrow(obs) - 1))
    flocker_data$ff_n_rep <- c(matrixStats::rowSums2(!is.na(obs)), 
                            rep(-99, nrow(obs) * (n_rep - 1)))
    flocker_data$ff_Q <- c(as.integer(matrixStats::rowSums2(obs, na.rm = T) > 0),
                        rep(-99, nrow(obs) * (n_rep - 1)))
    
    # Prepare to add rep indices, and trim flocker_data to existing observations
    flocker_data$ff_unit <- 1:nrow(obs)
    
    rep_indices <- as.data.frame(matrix(data = -99, nrow = nrow(flocker_data),
                                      ncol = n_rep))
    names(rep_indices) <- paste0("ff_rep_index", 1:n_rep)
    is_not_na <- !is.na(flocker_data$ff_y)
    rep_index_vec <- rep(-99, n_rep*nrow(obs))
    rep_index_vec[is_not_na] <- cumsum(is_not_na)[is_not_na]
    rep_indices[1:nrow(obs),] <- rep_index_vec
    
    flocker_data <- flocker_data[!is.na(flocker_data$ff_y), ]
    rep_indices <- rep_indices[is_not_na,]
    flocker_data <- cbind(flocker_data, rep_indices)
    
    out <- list(data = flocker_data, n_rep = n_rep,
                type = "single")
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
#' @return A flocker_data list that can be passed as data to \code{flock()}.
#' @export
make_flocker_data_dynamic <- function(obs, unit_covs = NULL, event_covs = NULL,
                                      quiet = FALSE) {
  n_year <- nslice(obs) # nslice checks that obs is a 3-D array
  n_series <- nrow(obs)
  n_rep <- ncol(obs)
  n_total <- n_year*n_series*n_rep
  
  assertthat::assert_that(
    n_year > 1, msg = "obs must contain at least two slices (seasons/years)"
  )
  assertthat::assert_that(
    n_rep > 1, 
    msg = paste0("obs must contain at least two columns (repeat visits to at ",
                 "least one unit)"
    )
  )
  unique_y <- unique(obs)[!is.na(unique(obs))]
  assertthat::assert_that(
    all(unique_y %in% c(0, 1)),
    msg = "obs may only contain the values 0, 1, or NA"
  )
  
  # Check that no NAs are non-trailing across columns (i.e. reps within series-
  # years)
  for (k in seq(n_year)) {
    for (j in 1:(n_rep - 1)) {
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
    !all(is.na(obs[ , n_rep, ])),
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
    assertthat::assert_that(
      is_named_list(event_covs),
      msg = "If provided, event_covs must be a named list."
    )
    n_event_covs <- length(event_covs)
    missing_covs <- vector()
    for (ec in 1:n_event_covs) {
      assertthat::assert_that(
        all.equal(dim(event_covs[[ec]]), dim(obs)),
        msg = paste0("Dimension mismatch found between obs and event_covs[[", ec, "]].")
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
  assertthat::assert_that(
    !is.null(event_covs),
    msg = paste0("Construction alert! The model contains no event covariates. ",
                 "This is fine, but for now please add a dummy event covariate.",
                 "You do not need to use this covariate in your model formula.")
  )
  
  # All unit covs are guaranteed to be not NA provided the unit is not part of, 
  # a block of trailing NAs, and all event covs are not NA provided that the 
  # response is not missing
  
  # add dummy data for the first visit to each unit, whether the unit was
  # sampled or not, unless the unit is part of a block of trailing NAs.
  n_year_obs <- apply(obs[ , 1, ], 1, max_position_not_na)
  assertthat::assert_that(
    !any(n_year_obs == 0),
    msg = paste0("at least one series (i.e. row; generally a site or a ",
                 "species-site) has no observations at any timestep")
  )
  for (i in seq_along(n_year_obs)) {
    obs[i, 1, 1 : n_year_obs[i]][is.na(obs[i, 1, 1 : n_year_obs[i]])] <- -99
  }
  
  # How many units will get modeled?
  rep1 <- obs[ , 1, ]
  n_unit <- sum(!is.na(rep1)) # We have already added a dummy value for rep1 if 
                              # it's NA and the unit is getting modeled.
  assertthat::assert_that(
    n_unit == sum(n_year_obs),
    msg = "error 4. This should not happen; please report a bug."
  )
  
  if (!is.null(event_covs)) {
    for (i in seq_along(event_covs)) {
      the_nas_event <- which(is.na(event_covs[[i]][ , 1, ]))
      the_m99s <- which(obs[ , 1, ] == -99)
      nas_replace <- shared_elements(the_nas_event, the_m99s)
      if (length(nas_replace) > 0) {
        if (!quiet) {
          warning(paste0("One or more event_covariates are NA on the first ",
                         "repeat visit to one or more units. This is ",
                         "allowable as long as the units in question were ",
                         "never visited (i.e. all visits NA in that site x ",
                         "timestep. However, to pass data to Stan, flocker will ",
                         "impute a covariate value corresponding to the mean of", 
                         "the event covariate across all sites/seasons.", 
                         "This imputation has no effect on the fitted posterior", 
                         "in any way, but could lead to spurious predicted", 
                         "detection probabilities on these never-conducted ", 
                         "visits."))
        }
        event_covs[[i]][ , 1, ][nas_replace] <- 
          mean(event_covs[[i]], na.rm = TRUE)
      }
    }
  }
  
  flocker_data <- data.frame(ff_y = expand_matrix(expand_array_3D(obs)))
  flocker_data$ff_n_unit <- c(n_unit, rep(-99, n_total - 1))
  if (!is.null(unit_covs)) {
    unit_covs_stacked <- stack_matrix(do.call(rbind, unit_covs), n_rep)
    flocker_data <- cbind(flocker_data, unit_covs_stacked)
  }
  event_covs <- as.data.frame(lapply(event_covs, function(x){expand_matrix(expand_array_3D((x)))}))
  flocker_data <- cbind(flocker_data, event_covs)
  # number of series (HMMs)
  flocker_data$ff_n_series <- c(n_series, rep(-99, n_total - 1))

  # For each series:
  # Number of units (seasons) in each series (don't marginalize over trailing
  # seasons with no observations)
  # all series have at least one unit
  flocker_data$ff_n_year <- c(apply(rep1, 1, max_position_not_na),
                           rep(-99, n_total - n_series))

  # For each unit:
  # Number of reps in each unit
  obs_unstacked <- expand_array_3D(obs)
  flocker_data$ff_n_rep <- c(
    apply(obs_unstacked, 1, max_position_not_na, treat_m99_NA = TRUE),
    rep(-99, n_total - n_series*n_year)
  )
  flocker_data$ff_Q <- c(as.integer(apply(obs_unstacked, 1, function(x){isTRUE(any(x == 1))})),
                      rep(-99, n_total - n_series*n_year))

  # Prepare to add unit and rep indices:
  flocker_data$ff_series <- seq(n_series)
  flocker_data$ff_year <- rep(seq(n_year), each = n_series)
  flocker_data$ff_series_year <- paste0(flocker_data$ff_series, "__", flocker_data$ff_year)
    
  # Trim flocker data to existing units
  flocker_data <- flocker_data[!is.na(flocker_data$ff_y), ]
  
  assertthat::assert_that(
    n_unit == length(unique(flocker_data$ff_series_year)),
    msg = "error 1; this should not happen; please report a bug"
  )
  
  # Unit indices: indices for the unique units belonging to each series
  unit_indices <- as.data.frame(matrix(data = -99, nrow = nrow(flocker_data),
                                      ncol = n_year))
  names(unit_indices) <- paste0("ff_unit_index", seq(n_year))
  if(!quiet){
    message("formatting unit indices")
    pb <- utils::txtProgressBar(min = 0, max = n_series, style = 3)
  }
  for (i in 1:n_series) {
    unit_indices[i, 1:flocker_data$ff_n_year[i]] <- 
      which(flocker_data$ff_series[seq(n_unit)] == flocker_data$ff_series[i])
    if(!quiet){
      utils::setTxtProgressBar(pb, i)
    }
  }
  if(!quiet){
    close(pb)
  }
  
  flocker_data <- cbind(flocker_data, unit_indices)
  
  # Rep indices: indices for the unique reps belonging to each unit
  rep_indices <- as.data.frame(matrix(data = -99, nrow = nrow(flocker_data),
                                       ncol = n_rep))
  names(rep_indices) <- paste0("ff_rep_index", seq(n_rep))
  
  if(!quiet){
    message("formatting rep indices")
    pb <- utils::txtProgressBar(min = 0, max = n_unit, style = 3)
  }

  for (i in seq(n_unit)) {
    assertthat::assert_that(
      !duplicated(flocker_data$ff_series_year)[i],
      msg = "error 2; this should not happen; please report a bug"
    )
    if (flocker_data$ff_n_rep[i] > 0) {
      rep_indices[i, 1:flocker_data$ff_n_rep[i]] <- 
        which(flocker_data$ff_series_year == flocker_data$ff_series_year[i])
    }
    if(!quiet){
      utils::setTxtProgressBar(pb, i)
    }
  }
  if(!quiet){
    close(pb)
  }
  
  flocker_data <- cbind(flocker_data, rep_indices)
  
  out <- list(data = flocker_data, n_rep = n_rep, n_year = n_year,
              type = "multi")
  class(out) <- c("list", "flocker_data")
  out
}

##### make_flocker_data_augmented #####

#' #' Format data for data-augmented occupancy model, to be passed to 
#'  \code{flock()}.
#' @param obs An I x J x K array where rows I are sites, columns J are 
#'  repeat sampling events, and slices K are species. Allowable values are 1 
#'  (detection), 0 (no detection), and NA (no sampling event).
#'   The data must be formatted so that all NAs are trailing within their rows.
#' @param n_aug Number of pseudospecies to augment
#' @param site_covs A dataframe of covariates for each site that are constant 
#'            across repeated sampling events.
#' @param event_covs A named list of I x J matrices, each one corresponding to a 
#' covariate that varies across repeated sampling events within sites
#' @param quiet Hide progress bars and informational messages?
#' @return A flocker_data list that can be passed as data to \code{flocker()}.
#' @export
make_flocker_data_augmented <- function(obs, n_aug, site_covs = NULL, 
                                        event_covs = NULL, quiet = FALSE) {
  assertthat::assert_that(
    length(dim(obs)) == 3,
    msg = "obs must have exactly three dimensions."
  )
  assertthat::assert_that(
    is_one_pos_int(n_aug),
    msg = "n_aug must be a positive integer"
  )
  obs1 <- obs[,,1]
  na_obs <- which(is.na(obs1))
  for (i in 2:dim(obs)[3]) {
    na_obs_i <- which(is.na(obs[,,i]))
    assertthat::assert_that(
      identical(na_obs, na_obs_i),
      msg = "Different species have different sampling events NA"
    )
  }
  n_rep <- ncol(obs1)
  assertthat::assert_that(
    n_rep >= 2,
    msg = "obs must contain at least two columns."
  )
  assertthat::assert_that(
    all(unique(obs) %in% c(0, 1, NA)),
    msg = "obs contains values other than 0, 1, NA."
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
      if (!identical(dim(event_covs[[ec]]), dim(obs1))) {
        stop(paste0("Dimension mismatch found between obs and event_covs[[", ec, "]]."))
      }
      missing_covs <- unique(c(missing_covs, which(is.na(event_covs[[ec]]))))
    }
    if (length(missing_covs) > 0) {
      if (!all(missing_covs %in% which(is.na(obs1)))) {
        stop(paste0("An event covariate contains missing values ",
                    "at a position where the response is not missing."))
      }
    }
  }
  
  n_sp_obs <- dim(obs)[3]
  n_sp <- n_sp_obs + n_aug
  n_site <- dim(obs)[1]
  aug_slice <- obs1
  aug_slice[!is.na(aug_slice)] <- 0
  dim(aug_slice) <- c(dim(aug_slice), 1)
  for (i in 1:n_aug) {
    obs <- abind::abind(obs, aug_slice, along = 3)
  }
  
  obs <- expand_array_3D(obs)
  
  flocker_data <- data.frame(ff_y = expand_matrix(obs))
  if (!is.null(site_covs)) {
    site_covs_stacked <- stack_matrix(site_covs, n_rep*n_sp)
    flocker_data <- cbind(flocker_data, site_covs_stacked)
  }
  if (!is.null(event_covs)) {
    event_covs <- lapply(event_covs, function(x){stack_matrix(x, n_sp)})
    event_covs <- as.data.frame(lapply(event_covs, expand_matrix))
    flocker_data <- cbind(flocker_data, event_covs)
  }
  
  flocker_data$ff_n_unit <- c(nrow(obs), 
                           rep(-99, nrow(obs) - 1))
  flocker_data$ff_n_rep <- c(apply(obs, 1, function(x){sum(!is.na(x))}), 
                          rep(-99, nrow(obs) * (n_rep - 1)))
  flocker_data$ff_Q <- c(as.integer(rowSums(obs, na.rm = T) > 0),
                      rep(-99, nrow(obs) * (n_rep - 1)))
  
  flocker_data$ff_n_sp <- c(n_sp, rep(-99, nrow(flocker_data)-1))
  flocker_data$ff_species <- rep(rep(c(1:n_sp), each = n_site), n_rep)
  flocker_data$ff_superQ <- c(rep(1, n_sp_obs), rep(0, n_aug), rep(-99, nrow(flocker_data) - n_sp))
  
  # Prepare to add rep indices, and trim flocker_data to existing observations
  flocker_data$ff_unit <- 1:nrow(obs)
  flocker_data <- flocker_data[!is.na(flocker_data$ff_y), ]
  rep_indices <- as.data.frame(matrix(data = -99, nrow = nrow(flocker_data),
                                      ncol = n_rep))
  names(rep_indices) <- paste0("ff_rep_index", 1:n_rep)
  if(!quiet){
    message("formatting rep indices")
    pb <- utils::txtProgressBar(min = 0, max = nrow(obs), style = 3)
  }
  for (i in 1:nrow(obs)) {
    rep_indices[i, 1:flocker_data$ff_n_rep[i]] <- which(flocker_data$ff_unit == i)
    if(!quiet){
      utils::setTxtProgressBar(pb, i)
    }
  }
  if(!quiet){
    close(pb)
  }
  flocker_data <- cbind(flocker_data, rep_indices)
  
  out <- list(data = flocker_data, n_rep = n_rep,
              type = "augmented")
  
  class(out) <- c("list", "flocker_data")
  out
}

