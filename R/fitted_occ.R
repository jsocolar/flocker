#' Get expected values of the posterior predictive distribution for the occupancy
#' component. Note that these are values for a hypothetical set of unit covariates
#' rather than the ones actually observed (some of which are known to be occupied,
#' when Q == 1).
#' @param flocker_fit A flocker_fit object.
#' @param new_data Optional new data at which to evaluate occupancy predictions. 
#'     If `NULL` (the default) expected values are generated for the original data.
#' @param CI A vector of length 2 specifying the upper and lower bounds of the 
#'     credible interval. Defaults to c(.05, .95)
#' @param ndraws Positive integer indicating how many poserior draws should be 
#'     used. If `NULL` (the default) all draws are used. 
#' @param response Should expected values be returned on the response or logit-scale? 
#' Defaults to `TRUE`, i.e. response scale.
#' @param re_formula formula containing group-level effects to be considered in 
#'     the prediction. If `NULL` (default), include all group-level effects; if 
#'     NA, include no group-level effects.

#' @return A data.frame containing the columns used to generate expected values 
#'     and the estimate and its credible interval (lower and upper bounds)
#' @export

fitted_occ <- function(flocker_fit, new_data = NULL, CI = c(.05, .95), 
                       ndraws = NULL, response=TRUE, re_formula = NULL, 
                       summarise = TRUE, allow_new_levels = FALSE) {
  ## catch errors
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
  
  linpred <- t(brms::posterior_linpred(flocker_fit, dpar = "occ", 
                                       ndraws = ndraws, 
                                       newdata = new_data_fmtd, 
                                       re_formula = re_formula))
  colnames(linpred) <- paste0("iter_", 1:ndraws)
  
  if(isTRUE(response)) {
    linpred <- boot::inv.logit(linpred)
  }
  
  if(isTRUE(summarise)) {
    linpred <- data.frame(estimate = matrixStats::rowMeans2(linpred),
                          lwr = matrixStats::rowQuantiles(linpred, probs = min(CI)), 
                          upr = matrixStats::rowQuantiles(linpred, probs = max(CI)))
    names(linpred)[2:3] <- paste0("Q", c(min(CI)*100, paste0(max(CI)*100)))
  }
  return(linpred)
}
