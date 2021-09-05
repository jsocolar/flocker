#' Fit an occupancy model
#' @param f_occ A brms-type model formula for occupancy. Must begin with "~".
#' @param f_det A brms-type model formula for detection. Must begin with "~".
#' @param flocker_data data, generally the output of \code{make_flocker_data()}.
#' @param data2 additional data (e.g. a covariance matrix for a phylogenetic effect)
#' @param ... additional arguments passed to \code{brms::brm()}
#' @param rep_constant A logical indicator. Are detection probabilities constant 
#'          across repeated sampling events within closure units?
#' @return a \code{brmsfit} containing the fitted occupancy model. 
#' @examples
#' \dontrun{
#' example_data <- example_flocker_data()
#' fd <- make_flocker_data(example_data$obs, example_data$unit_covs, example_data$event_covs)
#' flock(f_occ = ~ uc1 + s(uc2) + (1|grp),
#'           f_det = ~ uc1 + ec1 + s(ec2) + (1|grp),
#'           flocker_data = fd,
#'           refresh = 50, chains = 1, warmup = 5, iter = 200,
#'           adapt_engaged = F, step_size = .05, max_treedepth = 5, seed = 123)
#' }
#' @export

flock <- function(f_occ, f_det, flocker_data, data2 = NULL, 
                  rep_constant = FALSE, ...){
  if (!is.logical(rep_constant)) {
    stop("rep_constant must be logical.")
  }
  if (rep_constant) {
    if (flocker_data$type != "C") {
      stop(paste("flocker_data is not formatted for a model with rep-constant",
                 "detection probabilities"))
    }
  } else {
    if (flocker_data$type != "V") {
      stop(paste("flocker_data is not formatted for a model with rep-varying",
                 "detection probabilities"))
    }
    extra_args <- list(...)
    if ("threads" %in% names(extra_args)) {
      if (extra_args$threads > 1) {
        stop(paste("multithreading not allowed for a model with rep-specific",
                   "detection probabilities."))
      }
    }
  }
  if(!inherits(f_occ, "formula")) {
    stop(strwrap("Error in formula: does the occupancy formula have the correct 
                 syntax? e.g. ~ a + b"))
  }
  if(!inherits(f_det, "formula")) {
    stop(strwrap("Error in formula: does the detection formula have the correct 
                 syntax? e.g. ~ c + d"))
  }
  
  f_occ_txt <- paste0(deparse(f_occ), collapse = "")
  f_det_txt <- paste0(deparse(f_det), collapse = "")
  
  if (flocker_data$type == "V") {
    max_rep <- flocker_data$max_rep
    vint_text <- paste0("rep_index", 1:max_rep, 
                        collapse = ", ")
    
    f_occ_use <- stats::as.formula(paste0("occ ", f_occ_txt))
    f_det_use <- stats::as.formula(
      paste0("y | vint(n_unit, n_rep, Q, ",
             vint_text, ") ", f_det_txt))
    f_use <- brms::bf(f_det_use, f_occ_use)
    
    stanvars <- brms::stanvar(scode = make_occupancy_V_lpmf(max_rep = max_rep), block = "functions")
    flocker_fit <- brms::brm(f_use, 
                             data = flocker_data$data,
                             family = occupancy_V(max_rep), 
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
