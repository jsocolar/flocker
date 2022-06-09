#' Fit an occupancy model
#' @param f_occ A brms-type model formula for occupancy. If provided, must begin 
#'  with "~".
#' @param f_det A brms-type model formula for detection. Must begin with "~".
#' @param flocker_data data, generally the output of \code{make_flocker_data()}.
#' @param data2 additional data (e.g. a covariance matrix for a phylogenetic 
#'  effect)
#' @param multiseason Must be NULL (the default) or one of "colex" or 
#'  "autologistic". If NULL, data must be formatted for a single-season model.
#'  Otherwise, the data must be formatted for a multiseason model. If "colex", a 
#'  colonization-extinction model will be fit, and `f_col` and `f_ex` must be 
#'  specified. If "autologistic", an autologistic model will be fit, and `f_col`
#'  and `f_ex` must both be NULL.
#' @param f_col A brms-type model formula for colonization in 
#'  colonization-extinction dynamic models. If provided, must begin with "~".
#' @param f_ex A brms-type model formula for extinction probabilities in 
#'  colonization-extinction dynamic models. If provided, must begin with "~".
#' @param colex_init Must be NULL unless the model is a colonization-extinction
#'  dynamic model, in which case must be either "explicit" or "equilibrium". 
#'  If "explicit", then `f_occ` must be provided to model 
#'  of year-1 occupancy probabilities. If"equilibrium", then `f_occ` must be set
#'  to `NULL` and the initial occupancy probabilities are assumed to 
#'  be the (possibly site-specific) equilibrium probabilities from the 
#'  colonization-extinction dynamics. In this case, time-varying covariates are 
#'  not allowed for colonization or extinction (but are still allowed for 
#'  detection).
#' @param f_auto A brms-type model formula for the autologistic offset 
#'  parameter (theta). If provided, must begin with "~".
#' @param augmented Logical. Must be TRUE if data are formatted for a 
#'  data-augmented multi-species model, and FALSE otherwise.
#' @param fp logical. If data are formatted for an fp model, must be set to
#'  TRUE, otherwise to FALSE.
#' @param ... additional arguments passed to \code{brms::brm()}
#' @return a \code{brmsfit} containing the fitted occupancy model. 
#' @examples
#' \dontrun{
#' example_data <- example_flocker_data()
#' fd <- make_flocker_data(example_data$obs, example_data$unit_covs, 
#'                         example_data$event_covs)
#' flock(f_occ = ~ uc1 + s(uc2) + (1|grp),
#'           f_det = ~ uc1 + ec1 + s(ec2) + (1|grp),
#'           flocker_data = fd,
#'           refresh = 50, chains = 1, warmup = 5, iter = 200,
#'           adapt_engaged = F, step_size = .05, max_treedepth = 5, seed = 123)
#' }
#' @export
flock <- function(f_occ, f_det, flocker_data, data2 = NULL, 
                   multiseason = NULL, f_col = NULL, f_ex = NULL, colex_init = NULL, f_auto = NULL,
                   augmented = FALSE, fp = FALSE,
                   ...) {
  flock_(stancode = FALSE, f_occ = f_occ, f_det = f_det, 
         flocker_data = flocker_data, data2 = data2, 
         multiseason = multiseason, f_col = f_col, f_ex = f_ex, 
         colex_init = colex_init, f_auto = f_auto,
         augmented = augmented, fp = fp, ... = ...)
}



#' Generate stan code for an occupancy model
#' @inheritParams flock
#' @export
flocker_stancode <- function(f_occ, f_det, flocker_data, data2 = NULL, 
                  multiseason = NULL, f_col = NULL, f_ex = NULL, colex_init = NULL, f_auto = NULL,
                  augmented = FALSE, fp = FALSE,
                  ...) {
  flock_(stancode = TRUE, f_occ = f_occ, f_det = f_det, 
         flocker_data = flocker_data, data2 = data2, 
         multiseason = multiseason, f_col = f_col, f_ex = f_ex, 
         colex_init = colex_init, f_auto = f_auto,
         augmented = augmented, fp = fp, ... = ...)
}

#' Fit an occupancy model or generate Stan code for a model 
#' @param stancode logical: output Stan code or fitted model (T) or fit model (F)
#' @param f_occ A brms-type model formula for occupancy. If provided, must begin 
#'  with "~".
#' @param f_det A brms-type model formula for detection. Must begin with "~".
#' @param flocker_data data, generally the output of \code{make_flocker_data()}.
#' @param data2 additional data (e.g. a covariance matrix for a phylogenetic 
#'  effect)
#' @param multiseason Must be NULL (the default) or one of "colex" or 
#'  "autologistic". If NULL, data must be formatted for a single-season model.
#'  Otherwise, the data must be formatted for a multiseason model. If "colex", a 
#'  colonization-extinction model will be fit, and `f_col` and `f_ex` must be 
#'  specified. If "autologistic", an autologistic model will be fit, and `f_col`
#'  and `f_ex` must both be NULL.
#' @param f_col A brms-type model formula for colonization in 
#'  colonization-extinction dynamic models. If provided, must begin with "~".
#' @param f_ex A brms-type model formula for extinction probabilities in 
#'  colonization-extinction dynamic models. If provided, must begin with "~".
#' @param colex_init Must be NULL unless the model is a colonization-extinction
#'  dynamic model, in which case must be either "explicit" or "equilibrium". 
#'  If "explicit", then `f_occ` must be provided to model 
#'  of year-1 occupancy probabilities. If"equilibrium", then `f_occ` must be set
#'  to `NULL` and the initial occupancy probabilities are assumed to 
#'  be the (possibly site-specific) equilibrium probabilities from the 
#'  colonization-extinction dynamics. In this case, time-varying covariates are 
#'  not allowed for colonization or extinction (but are still allowed for 
#'  detection).
#' @param f_auto A brms-type model formula for the autologistic offset 
#'  parameter (theta). If provided, must begin with "~".
#' @param augmented Logical. Must be TRUE if data are formatted for a 
#'  data-augmented multi-species model, and FALSE otherwise.
#' @param fp logical. If data are formatted for an fp model, must be set to
#'  TRUE, otherwise to FALSE.
#' @param ... additional arguments passed to \code{brms::brm()}
#' @return stan code for brms model or a \code{brmsfit} containing the fitted occupancy model. 
flock_ <- function(stancode, f_occ, f_det, flocker_data, data2 = NULL, 
                  multiseason = NULL, f_col = NULL, f_ex = NULL, colex_init = NULL, f_auto = NULL,
                  augmented = FALSE, fp = FALSE,
                  ...) {
  ### validate parameters
  validate_flock_params(f_occ, f_det, flocker_data, multiseason, f_col, 
                        f_ex, colex_init, f_auto, augmented, fp)
  extra_args <- list(...)
  if ("threads" %in% names(extra_args)) {
    assertthat::assert_that(flocker_data$type %in% threading_types(),
    msg = "multithreading is not enabled for this model type")
  }
  
  ### create final formulas
  if (!is.null(f_occ)) {
    f_occ_txt <- paste0(deparse(f_occ), collapse = "")
    f_occ_use <- stats::as.formula(paste0("occ ", f_occ_txt))
  }
  if (!is.null(f_col)) {
    f_col_txt <- paste0(deparse(f_col), collapse = "")
    f_col_use <- stats::as.formula(paste0("colo ", f_col_txt))
  }
  if (!is.null(f_ex)) {
    f_ex_txt <- paste0(deparse(f_ex), collapse = "")
    f_ex_use <- stats::as.formula(paste0("ex ", f_ex_txt))
  }
  if (!is.null(f_auto)) {
    f_auto_txt <- paste0(deparse(f_auto), collapse = "")
    f_auto_use <- stats::as.formula(paste0("autologistic ", f_auto_txt))
  }
  
  f_det_txt <- paste0(deparse(f_det), collapse = "")
  # f_det_use is type specific and is handled below
  
  if (flocker_data$type == "single" & !fp) {
    max_rep <- flocker_data$n_rep
    vint_text <- paste0("rep_index", 1:max_rep, 
                        collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("y | vint(n_unit, n_rep, Q, ",
             vint_text, ") ", f_det_txt))
    f_use <- brms::bf(f_det_use, f_occ_use)
    stanvars <- brms::stanvar(
      scode = make_occupancy_single_lpmf(max_rep = max_rep), block = "functions"
    )
    out <- flocker_fit_code_util(stancode, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_single(max_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (flocker_data$type == "single_C") {
    f_det_use <- stats::as.formula(paste0("n_suc | vint(n_trial) ", f_det_txt))
    f_use <- brms::bf(f_det_use, f_occ_use)
    stanvars <- brms::stanvar(
      scode = make_occupancy_single_C_lpmf(), block = "functions"
    )
    out <- flocker_fit_code_util(stancode, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_single_C(), 
                                 stanvars = stanvars,
                                 ...)
  } else if (augmented) {
    max_rep <- flocker_data$n_rep
    vint_text <- paste0("rep_index", 1:max_rep, 
                        collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("y | vint(n_unit, n_rep, Q, n_sp, superQ, species, ",
             vint_text, ") ", f_det_txt))
    f_Omega_use <- stats::as.formula("Omega ~ 1")
    f_use <- brms::bf(f_det_use, f_occ_use, f_Omega_use)
    stanvars <- brms::stanvar(
      scode = make_occupancy_augmented_lpmf(max_rep = max_rep), 
      block = "functions"
    )
    out <- flocker_fit_code_util(stancode, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_augmented(max_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (isTRUE(multiseason == "colex") & isTRUE(colex_init == "explicit") & !fp) {
    n_rep <- flocker_data$n_rep
    n_year <- flocker_data$n_year
    vint_text1 <- paste0("unit_index", seq(n_year), 
                         collapse = ", ")
    vint_text2 <- paste0("rep_index", seq(n_rep), 
                        collapse = ", ")
    
    f_det_use <- stats::as.formula(
      paste0("y | vint(n_series, n_unit, n_year, n_rep, Q, ",
             vint_text1, ", ", vint_text2, ") ", f_det_txt))
    
    f_use <- brms::bf(f_det_use, f_occ_use, f_col_use, f_ex_use)
    
    stanvars <- 
      brms::stanvar(
        scode = make_emission_1(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_colex_likelihoods(), 
        block = "functions"
        ) + 
      brms::stanvar(
        scode = make_forward_colex(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_occupancy_multi_colex_lpmf(max_rep = n_rep, max_year = n_year), 
        block = "functions"
      )
    out <- flocker_fit_code_util(stancode, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_multi_colex(n_year, n_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (isTRUE(multiseason == "colex") & isTRUE(colex_init == "equilibrium") & !fp) {
    n_rep <- flocker_data$n_rep
    n_year <- flocker_data$n_year
    vint_text1 <- paste0("unit_index", seq(n_year), 
                         collapse = ", ")
    vint_text2 <- paste0("rep_index", seq(n_rep), 
                         collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("y | vint(n_series, n_unit, n_year, n_rep, Q, ",
             vint_text1, ", ", vint_text2, ") ", f_det_txt))
    
    f_use <- brms::bf(f_det_use, f_col_use, f_ex_use)
    
    stanvars <- 
      brms::stanvar(
        scode = make_emission_1(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_colex_likelihoods(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_forward_colex(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_occupancy_multi_colex_eq_lpmf(max_rep = n_rep, max_year = n_year), 
        block = "functions"
      )
    out <- flocker_fit_code_util(stancode, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_multi_colex_eq(n_year, n_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (isTRUE(multiseason == "autologistic") & !fp) {
    n_rep <- flocker_data$n_rep
    n_year <- flocker_data$n_year
    vint_text1 <- paste0("unit_index", seq(n_year), 
                         collapse = ", ")
    vint_text2 <- paste0("rep_index", seq(n_rep), 
                         collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("y | vint(n_series, n_unit, n_year, n_rep, Q, ",
             vint_text1, ", ", vint_text2, ") ", f_det_txt))
    
    f_use <- brms::bf(f_det_use, f_occ_use, f_auto_use)
    
    stanvars <- 
      brms::stanvar(
        scode = make_emission_1(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_colex_likelihoods(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_forward_colex(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_occupancy_multi_autologistic_lpmf(max_rep = n_rep, max_year = n_year), 
        block = "functions"
      )
    out <- flocker_fit_code_util(stancode, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_multi_autologistic(n_year, n_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (flocker_data$type == "single" & fp) {
    max_rep <- flocker_data$n_rep
    vint_text <- paste0("rep_index", 1:max_rep, 
                        collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("y | vint(n_unit, n_rep, Q, ",
             vint_text, ") ", f_det_txt))
    f_use <- brms::bf(f_det_use, f_occ_use)
    
    stanvars <- 
      brms::stanvar(scode = make_emission_0_fp(), block = "functions") +
      brms::stanvar(scode = make_emission_1_fp(), block = "functions") +
      brms::stanvar(
        scode = make_occupancy_single_fp_lpdf(max_rep = max_rep), 
        block = "functions"
      )
    out <- flocker_fit_code_util(stancode, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_single_fp(max_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (isTRUE(multiseason == "colex") & isTRUE(colex_init == "explicit") & fp) {
    n_rep <- flocker_data$n_rep
    n_year <- flocker_data$n_year
    vint_text1 <- paste0("unit_index", seq(n_year), 
                         collapse = ", ")
    vint_text2 <- paste0("rep_index", seq(n_rep), 
                         collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("y | vint(n_series, n_unit, n_year, n_rep, Q, ",
             vint_text1, ", ", vint_text2, ") ", f_det_txt))
    
    f_use <- brms::bf(f_det_use, f_occ_use, f_col_use, f_ex_use)
    
    stanvars <- 
      brms::stanvar(
        scode = make_emission_1_fp(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_emission_0_fp(),
        block = "functions"
      ) +
      brms::stanvar(
        scode = make_colex_fp_likelihoods(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_forward_colex_fp(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_occupancy_multi_colex_fp_lpdf(max_rep = n_rep, max_year = n_year), 
        block = "functions"
      )
    out <- flocker_fit_code_util(stancode, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_multi_colex_fp(n_year, n_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (isTRUE(multiseason == "colex") & isTRUE(colex_init == "equilibrium") & fp) {
    n_rep <- flocker_data$n_rep
    n_year <- flocker_data$n_year
    vint_text1 <- paste0("unit_index", seq(n_year), 
                         collapse = ", ")
    vint_text2 <- paste0("rep_index", seq(n_rep), 
                         collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("y | vint(n_series, n_unit, n_year, n_rep, Q, ",
             vint_text1, ", ", vint_text2, ") ", f_det_txt))
    
    f_use <- brms::bf(f_det_use, f_col_use, f_ex_use)
    
    stanvars <- 
      brms::stanvar(
        scode = make_emission_1_fp(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_emission_0_fp(),
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_colex_fp_likelihoods(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_forward_colex_fp(), 
        block = "functions"
      ) + 
      brms::stanvar(
        scode = make_occupancy_multi_colex_eq_fp_lpdf(max_rep = n_rep, max_year = n_year), 
        block = "functions"
      )
    out <- flocker_fit_code_util(stancode, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_multi_colex_eq_fp(n_year, n_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (isTRUE(multiseason == "autologistic") & fp) {
    stop("autologistic fp models not yet implemented")
  } else {
    stop("Error: unhandled type. This should not happen; please report a bug.")
  }
  if (!stancode) {
    attr(out, "class") <- c(attr(out, "class"), "flocker_fit")
    attr(out, "lik_type") <- flocker_data$type
    attr(out, "colex_init") <- flocker_data$colex_init
    attr(out, "fp") <- fp
  }
  out
}

#' Either fit a flocker model or output stancode
#' @param stancode logical; if TRUE return stan code, if false fit model
#' @param f_use brms formula
#' @param data brms data
#' @param data2 brms data2
#' @param family brms family
#' @param stanvars brms stanvars
#' @param ... brms ...
#' @return output of brms::make_stancode or brms::brm
flocker_fit_code_util <- function (
  stancode, f_use, data, data2, family, stanvars, ...
 ) {
  if (stancode) {
    out <- brms::make_stancode(f_use, 
                               data = data,
                               data2 = data2,
                               family = family, 
                               stanvars = stanvars,
                               ...)
  } else {
    out <- brms::brm(f_use, 
                     data = data,
                     data2 = data2,
                     family = family, 
                     stanvars = stanvars,
                     ...)
  }
  out
}
