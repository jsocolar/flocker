#' Get expected values of the posterior predictive distribution for the occupancy
#' component, the detection component, or both. Note that in the case of the 
#' detection component, you are conditioning on occupancy being 1. Note that these 
#' are values for a hypothetical set of unit covariates rather than the ones 
#' actually observed (some of which are known to be occupied, when Q == 1).
#' @param flocker_fit A flocker_fit object.
#' @param type Get posterior probabilities for 'occupancy', 'detection', or for 
#' 'both'. 
#' @param new_data Optional new data at which to evaluate occupancy predictions. 
#'     If `NULL` (the default) expected values are generated for the original data.
#' @param CI A vector of length 2 specifying the upper and lower bounds of the 
#'     credible interval. Defaults to c(.05, .95)
#' @param ndraws Positive integer indicating how many poserior draws should be 
#'     used. If `NULL` (the default) all draws are used. 
#' @param response Should expected values be returned on the response or logit-scale? 
#' Defaults to `TRUE`, i.e. response scale.
#' @param summarise if TRUE, return the expected value and upper and lower bound 
#' of the credible interval, otherwise return the posterior matrix. 
#' @param re_formula formula containing group-level effects to be considered in 
#'     the prediction. If `NULL` (default), include all group-level effects; if 
#'     NA, include no group-level effects.
#' @param sample_new_levels If new_data is provided and contains random effect
#'     levels not present in the original data, how should predictions be
#'     handled? Passed directly to brms::prepare_predictions.
#' @return A data.frame containing the columns used to generate expected values 
#'     and the estimate and its credible interval (lower and upper bounds)
#' @export
#' 
fitted_flocker <- function(flocker_fit, type, new_data = NULL, CI = c(.05, .95), 
                         ndraws = NULL, response=TRUE, re_formula = NULL, 
                         summarise = TRUE, allow_new_levels = FALSE, 
                         sample_new_levels = "uncertainty") {
  ## catch errors
  if(!(type %in% c("occupancy", "detection", "both"))) {
    stop("type must either be 'occupancy', 'detection', or 'both'")
  }
  if(length(CI) != 2 | !is.numeric(CI)) {
    stop("CI should be numeric length 2")
  }
  if(any(CI < 0 | CI > 1)) {
    stop("CI cannot have bounds <0 or >1")
  }
  if(!is.null(ndraws)) {
    if(ndraws > brms::ndraws(flocker_fit)) {
      stop("You can't have more draws than there are samples")
    }
  } else { # i.e. is null
    ndraws <- brms::ndraws(flocker_fit)
  }
  
  # format new_data: add necessary flocker cols to data frame
  if(is.null(new_data)) { # new data not provided
    new_data_fmtd <- flocker_fit$data
  } else {
    # add cols (not relevant for getting linpred on occ, but throws error 
    # otherwise)
    col_string <- "^y$|^Q$|^n_unit$|^unit$|^n_rep$|^rep_index"
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
  
  if(type %in% c("occupancy", "both")) {
    linpred_occ <- t(brms::posterior_linpred(flocker_fit, dpar = "occ", 
                                             ndraws = ndraws, 
                                             newdata = new_data_fmtd, 
                                             re_formula = re_formula, 
                                             allow_new_levels = allow_new_levels,
                                             sample_new_levels = sample_new_levels)) 
    linpred_occ <- boot::inv.logit(linpred_occ)
    if(type == "occupancy") linpred <- linpred_occ
  }
  if(type %in% c("detection", "both")) {
    linpred_det <- t(brms::posterior_linpred(flocker_fit, dpar = "mu", 
                                             ndraws = ndraws, 
                                             newdata = new_data_fmtd, 
                                             re_formula = re_formula, 
                                             allow_new_levels = allow_new_levels,
                                             sample_new_levels = sample_new_levels))
    linpred_det <- boot::inv.logit(linpred_det)
    if(type == "detection") linpred <- linpred_det
  }
  
  if(type == "both") {
    linpred <- exp(log(linpred_occ) + log(linpred_det))
  }
  
  colnames(linpred) <- paste0("iter_", 1:ndraws)
  
  if(isFALSE(response)) {
    linpred <- boot::logit(linpred)
  }
  
  if(isTRUE(summarise)) {
    linpred <- data.frame(estimate = matrixStats::rowMeans2(linpred),
                          lwr = matrixStats::rowQuantiles(linpred, probs = min(CI)), 
                          upr = matrixStats::rowQuantiles(linpred, probs = max(CI)))
    names(linpred)[2:3] <- paste0("Q", c(min(CI)*100, paste0(max(CI)*100)))
  }
  return(linpred)
}
