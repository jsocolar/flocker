#' Fit an occupancy model
#' @param f_occ A brms-type model formula for occupancy. If provided, must begin 
#'  with "~".
#' @param f_det A brms-type model formula for detection. Must begin with "~". 
#'  OR, a \code{brmsformula} object including formulas for all of the relevant
#'  distributional parameters in the desired model (det and at least one of occ, 
#'  colo, ex, autologistic, and Omega). The \code{$formula} element of the 
#'  \code{brmsformula} must be the detection formula, beginning with \code{det ~}.
#'  This latter option unadvisable except when necessary (e.g. when a nonlinear 
#'  formula is desired), as input checking is less thorough.
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
#' @param multi_init Must be NULL unless the model is a dynamic (multiseason)
#'  model, in which case must be either "explicit" or "equilibrium". 
#'  If "explicit", then `f_occ` must be provided to model occupancy
#'  probabilities in the first timestep. If "equilibrium", then `f_occ` must be
#'  `NULL` and the initial occupancy probabilities are assumed to be the 
#'  (possibly site-specific) equilibrium probabilities from the colonization-
#'  extinction dynamics.
#' @param f_auto Relevant only for autologistic models. A brms-type model 
#'   formula for the autologistic offset parameter (theta). If provided, must 
#'   begin with "~".
#' @param augmented Logical. Must be TRUE if data are formatted for a 
#'  data-augmented multi-species model, and FALSE otherwise.
#' @param threads NULL or positive integer. If integer, the number of threads
#'  to use per chain in within chain parallelization. Currently available only
#'  with single-season rep-constant models, and must be set to NULL otherwise.
#' @param ... additional arguments passed to \code{brms::brm()}
#' @return a \code{brmsfit} containing the fitted occupancy model. 
#' @examples
#' \donttest{
#' sfd <- simulate_flocker_data()
#' fd <- make_flocker_data(
#'  sfd$obs, 
#'  sfd$unit_covs,
#'  sfd$event_covs
#' )
#' flock(
#'   f_occ = ~ s(uc1) + (1|species),
#'   f_det = ~ uc1 + ec1 + (1|species),
#'   flocker_data = fd,
#'   refresh = 50, chains = 1, warmup = 5, iter = 200,
#'   control = list(adapt_engaged = FALSE, stepsize = .05, max_treedepth = 5),
#'   seed = 123
#'   )
#' }
#' @export
flock <- function(f_occ = NULL, f_det, flocker_data, data2 = NULL, 
                  multiseason = NULL, f_col = NULL, f_ex = NULL, 
                  multi_init = NULL, f_auto = NULL,
                  augmented = FALSE, threads = NULL,
                  ...) {
  flock_(output = "model", f_occ = f_occ, f_det = f_det, 
         flocker_data = flocker_data, data2 = data2, 
         multiseason = multiseason, f_col = f_col, f_ex = f_ex, 
         multi_init = multi_init, f_auto = f_auto,
         augmented = augmented, threads = threads, ...)
}

#' Generate stan code for an occupancy model
#' @inheritParams flock
#' @return generated stancode
#' @export
#' @examples
#' \donttest{
#' sfd <- simulate_flocker_data()
#' fd <- make_flocker_data(
#'  sfd$obs, 
#'  sfd$unit_covs,
#'  sfd$event_covs
#' )
#' flocker_stancode(
#'   f_occ = ~ s(uc1) + (1|species),
#'   f_det = ~ uc1 + ec1 + (1|species),
#'   flocker_data = fd,
#'   refresh = 50, chains = 1, warmup = 5, iter = 200,
#'   control = list(adapt_engaged = FALSE, stepsize = .05, max_treedepth = 5),
#'   seed = 123
#'   )
#'   }
flocker_stancode <- function(f_occ = NULL, f_det, flocker_data, data2 = NULL, 
                  multiseason = NULL, f_col = NULL, f_ex = NULL, multi_init = NULL, f_auto = NULL,
                  augmented = FALSE, threads = NULL,
                  ...) {
  flock_(output = "code", f_occ = f_occ, f_det = f_det, 
         flocker_data = flocker_data, data2 = data2, 
         multiseason = multiseason, f_col = f_col, f_ex = f_ex, 
         multi_init = multi_init, f_auto = f_auto,
         augmented = augmented, threads = threads, ...)
}

#' Generate stan data for an occupancy model
#' @inheritParams flock
#' @return generated stan data
#' @export
#' @examples
#' sfd <- simulate_flocker_data()
#' fd <- make_flocker_data(
#'  sfd$obs, 
#'  sfd$unit_covs,
#'  sfd$event_covs
#' )
#' flocker_standata(
#'   f_occ = ~ s(uc1) + (1|species),
#'   f_det = ~ uc1 + ec1 + (1|species),
#'   flocker_data = fd,
#'   refresh = 50, chains = 1, warmup = 5, iter = 200,
#'   control = list(adapt_engaged = FALSE, stepsize = .05, max_treedepth = 5),
#'   seed = 123
#'   )
flocker_standata <- function(f_occ=NULL, f_det, flocker_data, data2 = NULL, 
                             multiseason = NULL, f_col = NULL, f_ex = NULL, multi_init = NULL, f_auto = NULL,
                             augmented = FALSE, threads = NULL,
                             ...) {
  flock_(output = "data", f_occ = f_occ, f_det = f_det, 
         flocker_data = flocker_data, data2 = data2, 
         multiseason = multiseason, f_col = f_col, f_ex = f_ex, 
         multi_init = multi_init, f_auto = f_auto,
         augmented = augmented, threads = threads, ...)
}

#' Get prior for occupancy model
#' @inheritParams flock
#' @return A dataframe summarizing the parameters on which priors can
#'  be specified and giving the default priors for those parameters.
#'  See \code{?brms::get_prior} for further details.
#' @export
#' @examples
#' sfd <- simulate_flocker_data()
#' fd <- make_flocker_data(
#'  sfd$obs, 
#'  sfd$unit_covs,
#'  sfd$event_covs
#' )
#' get_flocker_prior(
#'  f_occ = ~ s(uc1) + + (1|species),
#'  f_det = ~ uc1 + ec1 + (1|species),
#'  flocker_data = fd
#' )
get_flocker_prior <- function(f_occ=NULL, f_det, flocker_data, data2 = NULL, 
                             multiseason = NULL, f_col = NULL, f_ex = NULL, multi_init = NULL, f_auto = NULL,
                             augmented = FALSE, threads = NULL,
                             ...) {
  flock_(output = "prior", f_occ = f_occ, f_det = f_det, 
         flocker_data = flocker_data, data2 = data2, 
         multiseason = multiseason, f_col = f_col, f_ex = f_ex, 
         multi_init = multi_init, f_auto = f_auto,
         augmented = augmented, threads = threads, ...)
}

#' Fit an occupancy model, or generate Stan code/data for a model, or get the prior
#' for a model
#' @param output "model" for model fitting, "code" for stancode, "data" for standata,
#'   "prior" for get_prior output
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
#' @param multi_init Must be NULL unless the model is a dynamic (multiseason)
#'  model, in which case must be either "explicit" or "equilibrium". 
#'  If "explicit", then `f_occ` must be provided to model occupancy
#'  probabilities in the first timestep. If "equilibrium", then `f_occ` must be
#'  `NULL` and the initial occupancy probabilities are assumed to be the 
#'  (possibly site-specific) equilibrium probabilities from the colonization-
#'  extinction dynamics.
#' @param f_auto A brms-type model formula for the autologistic offset 
#'  parameter (theta). If provided, must begin with "~".
#' @param augmented Logical. Must be TRUE if data are formatted for a 
#'  data-augmented multi-species model, and FALSE otherwise.
#' @param threads NULL or positive integer. If integer, the number of threads
#'  to use per chain in within chain parallelization.
#' @param ... additional arguments passed to \code{brms::brm()}
#' @return stan code for brms model or a \code{brmsfit} containing the fitted occupancy model.
#' @noRd
flock_ <- function(output, f_occ, f_det, flocker_data, data2 = NULL, 
                  multiseason = NULL, f_col = NULL, f_ex = NULL, multi_init = NULL, f_auto = NULL,
                  augmented = FALSE, threads = NULL,
                  ...) {
  ### validate parameters
  validate_flock_params(f_occ, f_det, flocker_data, multiseason, f_col, 
                        f_ex, multi_init, f_auto, augmented, threads)
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
  
  if(brms::is.brmsformula(f_det)){
    f_det_txt1 <- format(f_det$formula)
    is_valid <- grepl("^det[[:blank:]]*~", f_det_txt1)
    assertthat::assert_that(
      is_one_logical(is_valid),
      msg = "Error in formula checking. This should not happen; please report a bug."
    )
    assertthat::assert_that(
      is_valid,
      msg = paste0("When f_det is a brmsformula, the $formula element ", 
                   "must be a detection formula beginning with `det ~`"
      )
    )
    f_det_txt <- paste0("~", strsplit(f_det_txt1, "~")[[1]][2])
  } else if (brms::is.mvbrmsformula(f_det)) {
    assertthat::assert_that("det" %in% names(f_det$forms))
    f_det_txt1 <- format(f_det$forms$det$formula)
    is_valid <- grepl("^det[[:blank:]]*~", f_det_txt1)
    assertthat::assert_that(
      is_one_logical(is_valid),
      msg = "Error in formula checking. This should not happen; please report a bug."
    )
    assertthat::assert_that(
      is_valid,
      msg = paste0("When f_det is a brmsformula, the $formula element ", 
                   "must be a detection formula beginning with `det ~`"
      )
    )
    f_det_txt <- paste0("~", strsplit(f_det_txt1, "~")[[1]][2])
  } else {
    f_det_txt <- paste0(deparse(f_det), collapse = "")
  }
  
  # f_det_use is type specific and is handled below
  
  if (flocker_data$type == "single") {
    max_rep <- flocker_data$n_rep
    vint_text <- paste0("ff_rep_index", 1:max_rep, 
                        collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("ff_y | vint(ff_n_unit, ff_n_rep, ff_Q, ",
             vint_text, ") ", f_det_txt))
    if(brms::is.brmsformula(f_det)){
      f_use <- f_det
      f_use$formula <- f_det_use
    } else if(brms::is.mvbrmsformula(f_det)){
      f_use <- f_det
      f_use$forms$det$formula <- f_det_use
    } else {
      f_use <- brms::bf(f_det_use, f_occ_use)
    }
    if (is.null(threads)) {
      stanvars <- brms::stanvar(
        scode = make_occupancy_single_lpmf(max_rep = max_rep), block = "functions"
      )
      out <- flocker_fit_code_util(output, 
                                   f_use, 
                                   data = flocker_data$data,
                                   data2 = data2,
                                   family = occupancy_single(max_rep), 
                                   stanvars = stanvars,
                                   ...)
    } else {
      stanvars <- 
        brms::stanvar(
          scode = make_occupancy_single_partial_sum(max_rep = max_rep), 
          block = "functions"
        ) +
        brms::stanvar(
          scode = make_occupancy_single_threaded_lpmf(
            max_rep = max_rep,
            grainsize = as.integer(ceiling(flocker_data$data$ff_n_unit[1] / threads))
          ), 
          block = "functions"
        )
      out <- flocker_fit_code_util(output, 
                                   f_use, 
                                   data = flocker_data$data,
                                   data2 = data2,
                                   family = occupancy_single_threaded(max_rep), 
                                   stanvars = stanvars,
                                   threads = threads,
                                   ...)
    }

  } else if (flocker_data$type == "single_C") {
    f_det_use <- stats::as.formula(paste0("ff_n_suc | vint(ff_n_trial) ", f_det_txt))
    if(brms::is.brmsformula(f_det)){
      f_use <- f_det
      f_use$formula <- f_det_use
    } else {
      f_use <- brms::bf(f_det_use, f_occ_use)
    }
    stanvars <- brms::stanvar(
      scode = make_occupancy_single_C_lpmf(), block = "functions"
    )
    out <- flocker_fit_code_util(output, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_single_C(), 
                                 stanvars = stanvars,
                                 threads = threads,
                                 ...)
  } else if (augmented) {
    max_rep <- flocker_data$n_rep
    vint_text <- paste0("ff_rep_index", 1:max_rep, 
                        collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("ff_y | vint(ff_n_unit, ff_n_rep, ff_Q, ff_n_sp, ff_superQ, ff_species, ",
             vint_text, ") ", f_det_txt))
    
    if(brms::is.brmsformula(f_det)){
      f_use <- f_det
      f_use$formula <- f_det_use
    } else {
      f_Omega_use <- stats::as.formula("Omega ~ 1")
      f_use <- brms::bf(f_det_use, f_occ_use, f_Omega_use)
    }
    stanvars <- brms::stanvar(
      scode = make_occupancy_augmented_lpmf(max_rep = max_rep), 
      block = "functions"
    )
    out <- flocker_fit_code_util(output, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_augmented(max_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (isTRUE(multiseason == "colex") & isTRUE(multi_init == "explicit")) {
    n_rep <- flocker_data$n_rep
    n_year <- flocker_data$n_year
    vint_text1 <- paste0("ff_unit_index", seq(n_year), 
                         collapse = ", ")
    vint_text2 <- paste0("ff_rep_index", seq(n_rep), 
                        collapse = ", ")
    
    f_det_use <- stats::as.formula(
      paste0("ff_y | vint(ff_n_series, ff_n_unit, ff_n_year, ff_n_rep, ff_Q, ",
             vint_text1, ", ", vint_text2, ") ", f_det_txt))
    
    if(brms::is.brmsformula(f_det)){
      f_use <- f_det
      f_use$formula <- f_det_use
    } else {
      f_use <- brms::bf(f_det_use, f_occ_use, f_col_use, f_ex_use)
    }
    
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
    out <- flocker_fit_code_util(output, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_multi_colex(n_year, n_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (isTRUE(multiseason == "colex") & isTRUE(multi_init == "equilibrium")) {
    n_rep <- flocker_data$n_rep
    n_year <- flocker_data$n_year
    vint_text1 <- paste0("ff_unit_index", seq(n_year), 
                         collapse = ", ")
    vint_text2 <- paste0("ff_rep_index", seq(n_rep), 
                         collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("ff_y | vint(ff_n_series, ff_n_unit, ff_n_year, ff_n_rep, ff_Q, ",
             vint_text1, ", ", vint_text2, ") ", f_det_txt))
    
    if(brms::is.brmsformula(f_det)){
      f_use <- f_det
      f_use$formula <- f_det_use
    } else {
      f_use <- brms::bf(f_det_use, f_col_use, f_ex_use)
    }
    
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
    out <- flocker_fit_code_util(output, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_multi_colex_eq(n_year, n_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (isTRUE(multiseason == "autologistic") & isTRUE(multi_init == "explicit")) {
    n_rep <- flocker_data$n_rep
    n_year <- flocker_data$n_year
    vint_text1 <- paste0("ff_unit_index", seq(n_year), 
                         collapse = ", ")
    vint_text2 <- paste0("ff_rep_index", seq(n_rep), 
                         collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("ff_y | vint(ff_n_series, ff_n_unit, ff_n_year, ff_n_rep, ff_Q, ",
             vint_text1, ", ", vint_text2, ") ", f_det_txt))
    
    if(brms::is.brmsformula(f_det)){
      f_use <- f_det
      f_use$formula <- f_det_use
    } else {
      f_use <- brms::bf(f_det_use, f_occ_use, f_col_use, f_auto_use)
    }

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
    out <- flocker_fit_code_util(output, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_multi_autologistic(n_year, n_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else if (isTRUE(multiseason == "autologistic") & isTRUE(multi_init == "equilibrium")) {
    n_rep <- flocker_data$n_rep
    n_year <- flocker_data$n_year
    vint_text1 <- paste0("ff_unit_index", seq(n_year), 
                         collapse = ", ")
    vint_text2 <- paste0("ff_rep_index", seq(n_rep), 
                         collapse = ", ")
    f_det_use <- stats::as.formula(
      paste0("ff_y | vint(ff_n_series, ff_n_unit, ff_n_year, ff_n_rep, ff_Q, ",
             vint_text1, ", ", vint_text2, ") ", f_det_txt))
    
    if(brms::is.brmsformula(f_det)){
      f_use <- f_det
      f_use$formula <- f_det_use
    } else {
      f_use <- brms::bf(f_det_use, f_col_use, f_auto_use)
    }
    
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
        scode = make_occupancy_multi_autologistic_eq_lpmf(max_rep = n_rep, max_year = n_year), 
        block = "functions"
      )
    out <- flocker_fit_code_util(output, 
                                 f_use, 
                                 data = flocker_data$data,
                                 data2 = data2,
                                 family = occupancy_multi_autologistic_eq(n_year, n_rep), 
                                 stanvars = stanvars,
                                 ...)
  } else {
    stop("Error: unhandled type. This should not happen; please report a bug.")
  }
  if (output == "model") {
    attr(out, "class") <- c(attr(out, "class"), "flocker_fit")
    attr(out, "data_type") <- flocker_data$type
    attr(out, "multiseason") <- multiseason
    attr(out, "multi_init") <- multi_init
  }
  out
}

#' Either fit a flocker model or output stancode or output standata
#' @param output one of "model", "code", "data", or "prior"
#' @param f_use brms formula
#' @param data brms data
#' @param data2 brms data2
#' @param family brms family
#' @param stanvars brms stanvars
#' @param threads brms threads
#' @param ... brms ...
#' @return output of brms::brm or brms::make_stancode or brms::make_standata
#' @noRd
flocker_fit_code_util <- function (
  output, f_use, data, data2, family, stanvars, threads = NULL, ...
 ) {
  if (!is.null(threads)){ # this covers the rep-constant case where
    # native brms threading is available
    if (output == "code") {
      out <- brms::make_stancode(f_use, 
                                 data = data,
                                 data2 = data2,
                                 family = family, 
                                 stanvars = stanvars,
                                 threads = threads,
                                 ...)
    } else if (output == "data") {
      out <- brms::make_standata(f_use, 
                                 data = data,
                                 data2 = data2,
                                 family = family, 
                                 stanvars = stanvars,
                                 threads = threads,
                                 ...)
    } else if (output == "model") {
      out <- brms::brm(
        f_use,
        data = data,
        data2 = data2,
        family = family,
        stanvars = stanvars,
        threads = threads,
        ...
      )
    } else if (output == "prior") {
      out <- brms::get_prior(
        f_use, 
        data = data,
        data2 = data2,
        family = family,
        stanvars = stanvars,
        threads = threads,
        ...
      )
    }
  } else {
    if (output == "code") {
      out <- brms::make_stancode(
        f_use,
        data = data,
        data2 = data2,
        family = family, 
        stanvars = stanvars,
        ...
      )
    } else if (output == "data") {
      out <- brms::make_standata(
        f_use,
        data = data,
        data2 = data2,
        family = family, 
        stanvars = stanvars,
        ...
      )
    } else if (output == "model") {
      out <- brms::brm(
        f_use,
        data = data,
        data2 = data2,
        family = family,
        stanvars = stanvars,
        ...
      )
    } else if (output == "prior") {
      out <- brms::get_prior(
        f_use,
        data = data,
        data2 = data2,
        family = family,
        stanvars = stanvars,
        ...
      )
    }
  }
  out
}
