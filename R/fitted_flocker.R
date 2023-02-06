#' Get expected values of the posterior predictive distribution for the modeled
#' probabilities (occupancy, detection, colonization, extinction, autologistic). 
#' Note that these are conditional probabilties (e.g. detection conditional on 
#' occupancy, colonization conditional on previous non-occupancy, etc).
#' These probabilities are not conditioned on the observed histories (e.g. the
#' occupancy probability is not fixed to one at sites with a detection; it is
#' estimated only based on the covariates)
#' @param flocker_fit A flocker_fit object.
#' @param components a character vector specifying one or more of "occ",
#'     "det", "col", "ex", or "auto" for which to obtain fitted values.
#' @param new_data Optional new data at which to evaluate occupancy predictions. 
#'     New data can be passed as a flocker_data object produced by 
#'     \code{make_flocker_data} or as a dataframe with one row per desired
#'     predicton. If `NULL` (the default) expected values are generated for the 
#'     original data as formatted by make_flocker_data.
#' @param summarise if TRUE, return the expected value and upper and lower bound 
#'     of the credible interval, otherwise return posterior draws. 
#' @param CI A vector of length 2 specifying the upper and lower bounds of the 
#'     credible interval.
#' @param ndraws Positive integer indicating how many posterior draws should be 
#'     used. If `NULL` (the default) all draws are used. 
#' @param response Should results be returned on the response or logit scale? 
#'     Defaults to `TRUE`, i.e. response scale. However, the autologistic
#'     parameter is not interpretable as a probability and is always returned
#'     on the logit scale regardless of the value of `response`
#' @param re_formula formula containing group-level effects to be considered in 
#'     the prediction. If `NULL` (default), include all group-level effects; if 
#'     NA, include no group-level effects.
#' @param allow_new_levels allow new levels for random effect terms in `new_data`?
#'     Will error if set to `FALSE` and new levels are provided in `new_data`.
#' @param sample_new_levels If new_data is provided and contains random effect
#'     levels not present in the original data, how should predictions be
#'     handled? Passed directly to `brms::prepare_predictions`, which see.
#' @param unit_level Logical; defaults to FALSE. Relevant only when `new_data`
#'     is not a dataframe (i.e. it is `NULL` or a flocker_data object), and useful
#'     only for multiseason models with missing seasons. If FALSE, returns in the 
#'     shape of the observation matrix/array with NAs for missing visits. If
#'     TRUE, returns in the shape of the first visit, and returns values for all
#'     units that are not part of a trailing block of never-visited units,
#'     including never-visited units that are part of series with subsequent 
#'     visits.
#' @return A list of sets of expected values (one per component). If `new_data` 
#'     is a dataframe, each element contains one row per row of `new_data`.
#'     Otherwise, returns in the shape of the observation matrix/array used 
#'     to format the flocker_data (but see `unit_level` parameter for further
#'     details).
#' @export
#' 
fitted_flocker <- function(
    flocker_fit, 
    components = c("occ", "det", "col", "ex", "auto"),
    new_data = NULL, unit_level = FALSE, 
    summarise = FALSE, CI = c(.05, .95), ndraws = NULL, 
    response=TRUE, re_formula = NULL, allow_new_levels = FALSE, 
    sample_new_levels = "uncertainty") {
  
  assertthat::assert_that(
    is_flocker_fit(flocker_fit),
    msg = "flocker_fit is corrupt or is not a flocker fit."
  )
  assertthat::assert_that(
    is_one_logical(unit_level),
    msg = "unit_level must be a single logical value"
  )
  assertthat::assert_that(
    is_one_logical(summarise),
    msg = "summarise must be a single logical value"
  )
  if(summarise){
    assertthat::assert_that(
      is.numeric(CI),
      msg = "CI must be numeric"
    )
    assertthat::assert_that(
      length(CI) == 2,
      msg = "CI must be of length 2"
    )
    assertthat::assert_that(
      all(CI >= 0) & all(CI <= 1),
      msg = "CI must be between zero and one inclusive"
    )
  }
  if(!is.null(ndraws)) {
    assertthat::assert_that(
      ndraws <= brms::ndraws(flocker_fit),
      msg = "More draws requested than contained in flocker_fit"
    )
  } else {
    ndraws <- brms::ndraws(flocker_fit)
  }
  
  if(!is.null(new_data)) {
    assertthat::assert_that(
      isTRUE(class(new_data == "data.frame")) | isTRUE(is_flocker_data(new_data)),
      msg = "new_data must be NULL, a dataframe, or a flocker_data object"
    )
    if(is_flocker_data(new_data)) {
      assertthat::assert_that(
        new_data$type == attributes(flocker_fit)$data_type,
        msg = "new_data is formatted for a different model type than flocker_fit"
      )
    }
  }
  
  model_type <- type_flocker_fit(flocker_fit)
  relevant_components <- params_by_type[[model_type]]
  if("col" %in% components){components[components == "col"] <- "colo"}
  
  use_components <- components[components %in% relevant_components]
  
  assertthat::assert_that(
    length(use_components) > 0,
    msg = "none of the requested components is relevant to the supplied model type"
  )
  
  old_data <- flocker_fit$data
  
  # format new_data: add necessary flocker cols to data frame
  if(is.null(new_data)) { # new data not provided
    new_data_fmtd <- old_data
  } else if(is_flocker_data(new_data)) {
    new_data_fmtd <- new_data$data
  } else {
    # add cols to avoid error
    col_string <- "^ff_y$|^ff_Q$|^ff_n_unit$|^ff_unit$|^ff_n_rep$|^ff_rep_index|^ff_n_series|^ff_series|^ff_n_year|^ff_year|^ff_series_year"
    flocker_cols <- flocker_fit$data[grepl(col_string, names(flocker_fit$data))]
    extra_cols <- which(!(names(flocker_cols) %in% names(new_data)))
    if(length(extra_cols) == 0) {
      new_data_fmtd <- new_data
    } else {
      new_data_fmtd <- cbind(new_data, 
                             do.call(rbind, replicate(nrow(new_data), flocker_fit$data[1, extra_cols], F)), 
                             row.names=NULL) 
    }
  }
  
  assertthat::assert_that(
    all(names(flocker_fit$data) %in% names(new_data_fmtd)),
    msg = "new_data is missing columns required by flocker_fit (some covariates are missing or mis-named)"
  )
  
  component_list <- list()
  
  if("occ" %in% use_components) {
    linpred_occ <- t(brms::posterior_linpred(flocker_fit, dpar = "occ", 
                                             ndraws = ndraws, 
                                             newdata = new_data_fmtd, 
                                             re_formula = re_formula, 
                                             allow_new_levels = allow_new_levels,
                                             sample_new_levels = sample_new_levels))
    if(response){
      linpred_occ <- boot::inv.logit(linpred_occ)
    }
    component_list$linpred_occ <- linpred_occ
  }
  if("colo" %in% use_components) {
    linpred_colo <- t(brms::posterior_linpred(flocker_fit, dpar = "colo", 
                                             ndraws = ndraws, 
                                             newdata = new_data_fmtd, 
                                             re_formula = re_formula, 
                                             allow_new_levels = allow_new_levels,
                                             sample_new_levels = sample_new_levels))
    if(response){
      linpred_colo <- boot::inv.logit(linpred_colo)
    }
    component_list$linpred_col <- linpred_colo
  }
  if("ex" %in% use_components) {
    linpred_ex <- t(brms::posterior_linpred(flocker_fit, dpar = "ex", 
                                             ndraws = ndraws, 
                                             newdata = new_data_fmtd, 
                                             re_formula = re_formula, 
                                             allow_new_levels = allow_new_levels,
                                             sample_new_levels = sample_new_levels))
    if(response){
      linpred_ex <- boot::inv.logit(linpred_ex)
    }
    component_list$linpred_ex <- linpred_ex
  }
  if("auto" %in% use_components) {
    linpred_auto <- t(brms::posterior_linpred(flocker_fit, dpar = "auto", 
                                             ndraws = ndraws, 
                                             newdata = new_data_fmtd, 
                                             re_formula = re_formula, 
                                             allow_new_levels = allow_new_levels,
                                             sample_new_levels = sample_new_levels))
    if(response){
      linpred_auto <- boot::inv.logit(linpred_auto)
    }
    component_list$linpred_auto <- linpred_auto
  }
  if("det" %in% use_components) {
    linpred_det <- t(brms::posterior_linpred(flocker_fit, dpar = "mu", 
                                             ndraws = ndraws, 
                                             newdata = new_data_fmtd, 
                                             re_formula = re_formula, 
                                             allow_new_levels = allow_new_levels,
                                             sample_new_levels = sample_new_levels))
    if(response){
      linpred_det <- boot::inv.logit(linpred_det)
    }
    component_list$linpred_det <- linpred_det
  }
  
  if(summarise) {
    cl2 <- lapply(component_list, summarise_fun)
  } else {
    cl2 <- component_list
  }
  
  if(is.null(new_data)) {
    gp <- get_positions(flocker_fit, unit_level = unit_level)
    out <- lapply(cl2, reshape_fun, gp = gp) ##HERE WE ARE
  } else if (is_flocker_data(new_data)) {
    gp <- get_positions(new_data, unit_level = unit_level)
    out <- lapply(cl2, reshape_fun, gp = gp) ##HERE WE ARE
  } else {
    out <- cl2
  }
    
  out
}

#' function to summarize matrix of linear predictors to its mean and CI
#' @param x linpreds
#' @return summary
summarise_fun <- function(x) {
  out <- data.frame(estimate = matrixStats::rowMeans2(x),
             lwr = matrixStats::rowQuantiles(x, probs = min(CI)), 
             upr = matrixStats::rowQuantiles(x, probs = max(CI)))
  names(out)[2:3] <- paste0("Q", c(min(CI)*100, paste0(max(CI)*100)))
  out
}

#' function to reshape a matrix of values into a stack of arrays via get_positions
#' @param x input_matrix
#' @param gp output of get_positions
reshape_fun <- function(x, gp) {
  assertthat::assert_that(is.matrix(x))
  ai <- list()
  for(i in seq_len(ncol(x))) {
    ai[[i]] <- array(x[,i][gp], dim = dim(gp))
  }
  if(ncol(x) > 1) {
    arr_dim <- length(dim(ai[[1]]))
    return(abind::abind(ai, along = arr_dim + 1))
  } else {
    return(ai[[1]]) 
  }
}
