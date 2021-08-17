#' Fit an occupancy model
#' @param f_occ A brms-type model formula for occupancy. Must begin with "~".
#'              # Important
#'              If \code{visit_constant = T}, then the occupancy sub-model gives the 
#'              probability of NON-occupancy. Negate the terms associated with this 
#'              sub-model to recover covariate effects on occupancy.
#' @param f_det A brms-type model formula for detection. Must begin with "~".
#' @param flocker_data data, generally the output of \code{make_flocker_data()}.
#' @param data2 additional data (e.g. a covariance matrix for a phylogenetic effect)
#' @param ... additional arguments passed to \code{brms::brm()}
#' @param visit_constant A logical indicator. Are detection probabilities constant across visits?
#'          Must be TRUE if the model lacks visit-specific detection covariates, otherwise must
#'          be FALSE.
#' @return a \code{brmsfit} containing the fitted occupancy model. 
#' @examples
#' \dontrun{
#' example_data <- example_flocker_data()
#' fd <- make_flocker_data(example_data$obs, example_data$unit_covs, example_data$visit_covs)
#' flock(f_occ = ~ sc1 + s(sc2) + (1|grp),
#'           f_det = ~ sc1 + vc1 + s(vc2) + (1|grp),
#'           flocker_data = fd,
#'           refresh = 50, chains = 1, warmup = 5, iter = 200,
#'           adapt_engaged = F, step_size = .05, max_treedepth = 5, seed = 123)
#' }
#' @export

flock <- function(f_occ, f_det, flocker_data, data2 = NULL, 
                  visit_constant = FALSE, ...){
  if (!is.logical(visit_constant)) {
    stop("visit_constant must be logical.")
  }
  if (visit_constant) {
    if (flocker_data$type != "C") {
      stop(paste("flocker_data is not formatted for a model with visit-constant",
                 "detection probabilities"))
    }
  } else {
    if (flocker_data$type != "V") {
      stop(paste("flocker_data is not formatted for a model with visit-specific",
                 "detection probabilities"))
    }
    extra_args <- list(...)
    if ("threads" %in% names(extra_args)) {
      if (extra_args$threads > 1) {
        stop(paste("multithreading not allowed for a model with visit-specific",
                    "detection probabilities."))
      }
    }
  }
  
  f_occ_txt <- tryCatch(paste0(deparse(f_occ), collapse = ""),
                        error = function(x) 
                          strwrap("Error in formula deparse: does the occupancy 
                                  formula have the correct syntax? e.g. ~ a + b"))
  
  f_det_txt <- tryCatch(paste0(deparse(f_det), collapse = ""), 
                        error = function(x) 
                          strwrap("Error in formula deparse: does the detection 
                                  formula have the correct syntax? e.g. ~ c + d"))
  
  if (flocker_data$type == "V") {
    max_visit <- flocker_data$max_visit
    vint_text <- paste0("visit_index", 1:max_visit, 
                        collapse = ", ")
    
    f_occ_use <- stats::as.formula(paste0("occ ", f_occ_txt))
    f_det_use <- stats::as.formula(
      paste0("y | vint(n_unit, n_visit, Q, ",
             vint_text, ") ", f_det_txt))
    f_use <- brms::bf(f_det_use, f_occ_use)
  
    stanvars <- brms::stanvar(scode = make_occupancy_V_lpmf(max_visit = max_visit), block = "functions")
    flocker_fit <- brms::brm(f_use, 
                             data = flocker_data$data,
                             family = occupancy_V(max_visit), 
                             stanvars = stanvars,
                             ...)
  } else if (flocker_data$type == "C") {
    f_occ_use <- stats::as.formula(paste0("occ ", f_occ_txt))
    f_det_use <- stats::as.formula(paste0("n_suc | vint(n_trial) ", f_det_txt))
    f_use <- brms::bf(f_det_use, f_occ_use)
    stanvars <- brms::stanvar(scode = make_occupancy_C_lpmf(), block = "functions")
    flocker_fit <- brms::brm(brms::bf(f_use),
                             data = flocker_data$data,
                             family = occupancy_C(),
                             stanvars = stanvars,
                             ...)
  }
  attr(flocker_fit, "class") <- c(attr(flocker_fit, "class"), "flocker_fit")
  attr(flocker_fit, "lik_type") <- flocker_data$type
  flocker_fit
}
