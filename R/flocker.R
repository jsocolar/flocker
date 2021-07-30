#' Fit an occupancy model
#' @param f_occ A brms-type model formula for occupancy. Must begin with "~".
#'              # Important
#'              If visit_constant = T, then the occupancy sub-model gives the probability
#'              of NON-occupancy. Negate the terms associated with this sub-model to recover 
#'              covariate effects on occupancy.
#' @param f_det A brms-type model formula for detection. Must begin with "~".
#' @param flocker_data data, generally the output of `make_flocker_data()`.
#' @param data2 additional data (e.g. a covariance matrix for a phylogenetic effect)
#' @param ... additional arguments passed to cmdstanr::sample() if visit_constant is FALSE
#'            or to brms::brm() if visit_constant is TRUE
#' @param visit_constant A logical indicator. Are detection probabilities constant across visits?
#'          Must be TRUE if the model lacks visit-specific detection covariates, otherwise must
#'          be FALSE.
#' @return the fitted occupancy model. If visit_constant = F, a cmdstan_fit object from cmdstanr.
#'         If visit_constant = T, a brmsfit object from brms.
#'         # Important
#'         If visit_constant = T, then the occupancy sub-model gives the probability
#'         of NON-occupancy. Negate the terms associated with this sub-model to recover covariate
#'         effects on occupancy.
#' @example
#' example_data <- example_flocker_data()
#' fd <- make_flocker_data(example_data$obs, example_data$site_covs, example_data$visit_covs)
#' flocker(f_occ = ~ sc1 + s(sc2) + (1|grp), 
#'         f_det = ~ sc1 + vc1 + s(vc2) + (1|grp), 
#'         flocker_data = fd, 
#'         refresh = 1, chains = 1, iter_warmup = 5, iter_sampling = 5)
#' @import cmdstanr
#' @export

flocker <- function(f_occ, f_det, flocker_data, data2 = NULL, visit_constant = FALSE, ...){
  if (!is.logical(visit_constant)) {
    stop("visit_constant must be logical")
  }
  if (visit_constant & (flocker_data$.type != "N")) {
    stop("flocker_data is not formatted for a model with visit-constant detection probabilities")
  }else if ((!visit_constant) & (flocker_data$.type != "V")) {
    stop("flocker_data is not formatted for a model with visit-specific detection probabilities")
  }
  
  f_occ_txt <- paste0(deparse(f_occ), collapse = "")
  f_det_txt <- paste0(deparse(f_det), collapse = "")
  
  if (flocker_data$.type == "V") {
    f_occ_use <- as.formula(paste0("occ | resp_subset(occupancy_subset) ", f_occ_txt))
    f_det_use <- as.formula(paste0("det ", f_det_txt))
    
    # make code and data
    flocker_stancode <- flocker_make_stancode(f_occ_use, f_det_use, flocker_data, data2)
    flocker_standata <- flocker_make_standata(f_occ_use, f_det_use, flocker_data, data2)
    
    # write .stan file in temp directory
    fileConn<-file(paste0(tempdir(), "/flocker_model.stan"))
    writeLines(flocker_stancode, fileConn)
    close(fileConn)
    
    # compile model
    flocker_model <- cmdstanr::cmdstan_model(paste0(tempdir(), "/flocker_model.stan"))
    
    # sample model
    flocker_fit <- flocker_model$sample(data = flocker_standata, ...)
  } else if (flocker_data$.type == "N") {
    f_occ_use <- as.formula(paste0("n_suc | trials(n_trial) ", f_occ_txt))
    f_det_use <- as.formula(paste0("zi ", f_det_txt))
    
    flocker_fit <- brms::brm(brms::bf(f_occ_use, f_det_use),
                             data = flocker_data$flocker_data,
                             family = brms::zero_inflated_binomial(),
                             backend = 'cmdstanr', ...)
  }
  
  flocker_fit
}
