#' Posterior predictive distributions for modeled probabilities
#' 
#' Get expected values of the posterior predictive distribution for the modeled
#' probabilities (occupancy, detection, colonization, extinction, autologistic).
#' 
#' The probabilities returned are conditional probabilities (e.g. detection 
#' conditional on occupancy, colonization conditional on previous 
#' non-occupancy, etc). These probabilities are not conditioned on the 
#' observed histories (e.g. the occupancy probability is not fixed to one 
#' at sites with a detection; it is estimated only based on the covariates).
#' @param flocker_fit A flocker_fit object.
#' @param components a character vector specifying one or more of "occ",
#'     "det", "col", "ex", "auto", and "Omega" for which to obtain fitted values.
#' @param new_data Optional new data at which to evaluate occupancy predictions. 
#'     New data can be passed as a flocker_data object produced by 
#'     \code{make_flocker_data} or as a dataframe with one row per desired
#'     prediction. If `NULL` (the default) expected values are generated for the 
#'     original data as formatted by make_flocker_data.
#' @param summarise if TRUE, return the expected value and upper and lower bound 
#'     of the credible interval, otherwise return posterior draws. 
#' @param CI A vector of length 2 specifying the upper and lower bounds of the 
#'     credible interval.
#' @param draw_ids Vector of indices of the posterior draws to be 
#'     used. If `NULL` (the default) all draws are used in their native order. 
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
#' @examples 
#' \donttest{
#' fitted_flocker(
#'   example_flocker_model_single, 
#'   summarise = TRUE
#' )
#' }
fitted_flocker <- function(
    flocker_fit, 
    components = c("occ", "det", "col", "ex", "auto", "Omega"),
    new_data = NULL, unit_level = FALSE, 
    summarise = FALSE, CI = c(.05, .95), draw_ids = NULL, 
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
  ndraws_total <- brms::ndraws(flocker_fit)
  if(!is.null(draw_ids)) {
    assertthat::assert_that(
      max(draw_ids) <= ndraws_total,
      msg = "A requested draw id is greater than the number of draws contained in flocker_fit"
    )
  } else {
    draw_ids <- seq_len(brms::ndraws(flocker_fit))
  }
  
  if(!is.null(new_data)) {
    assertthat::assert_that(
      isTRUE(class(new_data) == "data.frame") | isTRUE(is_flocker_data(new_data)),
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
  
  if(is.null(new_data)) { # new data not provided
    new_data_fmtd <- old_data
  } else if(is_flocker_data(new_data)) {
    new_data_fmtd <- new_data$data
  } else {
    new_data_fmtd <- new_data
  }
  
  # Handle models fitted via 0 + Intercept syntax with newdata
  if(
    ("Intercept" %in% names(flocker_fit$data)) &
    !("Intercept" %in% names(new_data_fmtd))){
    new_data_fmtd$Intercept <- 1
  }
  
  f_form <- flocker_fit$formula
  
  component_list <- list()
  
  if("occ" %in% use_components) {
    occ_covs <- all.vars(f_form[[2]]$occ[[3]])
    linpred_occ <- t(
      brms::posterior_linpred(
        flocker_fit, dpar = "occ", draw_ids = draw_ids, newdata = new_data_fmtd, 
        re_formula = re_formula, allow_new_levels = allow_new_levels, 
        sample_new_levels = sample_new_levels,
        req_vars = occ_covs
        )
      )
    if(response){
      linpred_occ <- boot::inv.logit(linpred_occ)
    }
    component_list$linpred_occ <- linpred_occ
  }
  if("colo" %in% use_components) {
    colo_covs <- all.vars(f_form[[2]]$colo[[3]])
    linpred_colo <- t(
      brms::posterior_linpred(
        flocker_fit, dpar = "colo", draw_ids = draw_ids, newdata = new_data_fmtd, 
        re_formula = re_formula, allow_new_levels = allow_new_levels, 
        sample_new_levels = sample_new_levels,
        req_vars = colo_covs
        )
      )
    if(response){
      linpred_colo <- boot::inv.logit(linpred_colo)
    }
    component_list$linpred_col <- linpred_colo
  }
  if("ex" %in% use_components) {
    ex_covs <- all.vars(f_form[[2]]$ex[[3]])
    linpred_ex <- t(
      brms::posterior_linpred(
        flocker_fit, dpar = "ex", draw_ids = draw_ids, newdata = new_data_fmtd, 
        re_formula = re_formula, allow_new_levels = allow_new_levels, 
        sample_new_levels = sample_new_levels,
        req_vars = ex_covs
        )
      )
    if(response){
      linpred_ex <- boot::inv.logit(linpred_ex)
    }
    component_list$linpred_ex <- linpred_ex
  }
  if("auto" %in% use_components) {
    auto_covs <- all.vars(f_form[[2]]$autologistic[[3]])
    linpred_auto <- t(
      brms::posterior_linpred(
        flocker_fit, dpar = "autologistic", draw_ids = draw_ids, newdata = new_data_fmtd, 
        re_formula = re_formula, allow_new_levels = allow_new_levels, 
        sample_new_levels = sample_new_levels,
        req_vars = auto_covs
        )
      )

    component_list$linpred_auto <- linpred_auto
  }
  if("det" %in% use_components) {
    det_covs <- all.vars(f_form[[1]][[3]])
    linpred_det <- t(
      brms::posterior_linpred(
        flocker_fit, dpar = "mu", draw_ids = draw_ids, newdata = new_data_fmtd, 
        re_formula = re_formula, allow_new_levels = allow_new_levels,
        sample_new_levels = sample_new_levels,
        req_vars = det_covs
        )
      )
    if(response){
      linpred_det <- boot::inv.logit(linpred_det)
    }
    component_list$linpred_det <- linpred_det
  }
  
  if("Omega" %in% use_components) {
    Omega_covs <- all.vars(f_form[[2]]$Omega[[3]])
    
    linpred_Omega <- t(
      brms::posterior_linpred(
        flocker_fit, dpar = "Omega", draw_ids = draw_ids, newdata = new_data_fmtd, 
        re_formula = re_formula, allow_new_levels = allow_new_levels, 
        sample_new_levels = sample_new_levels,
        req_vars = Omega_covs
        )
      )
    if(response){
      linpred_Omega <- boot::inv.logit(linpred_Omega)
    }
    component_list$linpred_Omega <- linpred_Omega
  }
  
  if(summarise) {
    cl2 <- lapply(component_list, summarise_fun, CI = CI)
  } else {
    cl2 <- component_list
  }
  
  if(is.null(new_data)) {
    gp <- get_positions(flocker_fit, unit_level = unit_level)
    out <- lapply(cl2, reshape_fun, gp = gp)
  } else if (is_flocker_data(new_data)) {
    gp <- get_positions(new_data, unit_level = unit_level)
    out <- lapply(cl2, reshape_fun, gp = gp)
  } else {
    out <- cl2
  }
  
  dt <- attributes(flocker_fit)$data_type
  
  if(is.null(new_data) | is_flocker_data(new_data)) {
    dn <- list(
      paste0("site_", seq_len(dim(out[[1]])[1]))
    )
    if(!unit_level) {
      dn <- append(dn, list(paste0("visit_", seq_len(dim(out[[1]])[2]))))
      multi_ind <- 3
    } else {
      multi_ind <- 2
    }
    if(dt == "augmented"){
      dn <- append(dn, list(paste0("species_", seq_len(dim(out[[1]])[multi_ind]))))
    }
    if(dt == "multi"){
      dn <- append(dn, list(paste0("timestep_", seq_len(dim(out[[1]])[multi_ind]))))
    }
  } else {
    dn <- list(
      paste0("row_", seq_len(dim(out[[1]])[1]))
    )
  }
  
  if(summarise){
    dn <- append(dn, list(c("mean", paste0("Q", c(min(CI)*100, paste0(max(CI)*100))))))
  } else {
    dn <- append(dn, list(paste0("draw_", draw_ids)))
  }
  
  for(i in seq_along(out)){
    dimnames(out[[i]]) <- dn
  }
  out
}

#' summarize matrix of linear predictors to its mean and CI
#' @param x linpreds
#' @param CI fractional bounds of credible interval between 0 and 1 inclusive
#' @return summary
#' @noRd
summarise_fun <- function(x, CI) {
  out <- data.frame(estimate = matrixStats::rowMeans2(x),
             lwr = matrixStats::rowQuantiles(x, probs = min(CI)), 
             upr = matrixStats::rowQuantiles(x, probs = max(CI)))
  names(out)[2:3] <- paste0("Q", c(min(CI)*100, paste0(max(CI)*100)))
  out
}

#' reshape a matrix of values into a stack of arrays via get_positions
#' @param x input_matrix
#' @param gp output of get_positions
#' @noRd
reshape_fun <- function(x, gp) {
  assertthat::assert_that(is.matrix(x) | is.data.frame(x))
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
