#' Get expected values of the posterior predictive distribution for the occupancy
#' component. Note that these are values for a hypothetical set of unit covariates
#' rather than the ones actually observed (some of which are known to be occupied,
#' when Q == 1).
#' @param flocker_fit A flocker_fit object.
#' @param new_data Optional new data at which to evaluate occupancy predictions. 
#'     If `NULL` (the default) expected values are generated for the original data.
#' @param CI A vector of length 2 specifying the upper and lower bounds of the 
#'     credible interval
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

# Note: if user doesn't specify new_data, the returned dataframe will contain 
# a number of irrelevant columns. Need to fix this. Also needs editing to work 
# for phylogenetic models (or anything using data2). 
fitted_occ <- function(flocker_fit, new_data = NULL, CI = c(.05, .95), ndraws = NULL, 
                       response=TRUE, re_formula = NULL) {
    # catch errors
    if(length(CI) != 2) {
        stop("CI should just be length 2 (lwr, upr)")
    }
    
    if(any(CI < 0 | CI > 1)) {
        stop("CI cannot have bounds <0 or >1")
    }
    
    if(!is.null(ndraws) & ndraws > brms::ndraws(flocker_fit)) {
        stop("You can't have more draws than there are samples")
    }
    if(is.null(ndraws)) {
        ndraws <- brms::ndraws(flocker_fit)
    }
    
    # format new_data
    if(is.null(new_data)) { 
        # only need first block for occ linpred
        new_data_fmtd <- flocker_fit$data[1:flocker_fit$data$n_unit[1],]
    } else {
        # add cols (not relevant for getting linpred on occ, but throws error 
        # otherwise)
        extra_cols <- which(!(names(fit_ed$data) %in% names(new_data)))
        new_data_fmtd <- cbind(new_data, flocker_fit$data[1, extra_cols], row.names=NULL)
    } 
    
    linpred <- t(brms::posterior_linpred(flocker_fit, dpar = "occ", 
                                         ndraws = ndraws, 
                                         newdata = new_data_fmtd, 
                                         re_formula = re_formula))
    if(response == TRUE) {
        linpred <- boot::inv.logit(linpred)
    }
    ests <- data.frame(estimate = matrixStats::rowMeans2(linpred),
                       lwr = matrixStats::rowQuantiles(linpred, probs = min(CI)), 
                       upr = matrixStats::rowQuantiles(linpred, probs = max(CI)))
    # return
    cbind(new_data, ests)
}