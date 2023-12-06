#' A posterior predict function for the rep-constant occupancy model, sufficient 
#' for \code{brms::predict.brmsfit()} to work. 
#' @param i Observation id
#' @param prep Output of \code{brms::prepare_predictions}. See brms custom 
#' families vignette at 
#' https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html
#' @param ... unused additional arguments. See brms custom 
#' families vignette at 
#' https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html
#' @return Posterior predictions
#' @noRd
posterior_predict_occupancy_single_C <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  occ <- brms::get_dpar(prep, "occ", i = i)
  trials <- prep$data$vint1[i]
  stats::rbinom(length(mu), 1, boot::inv.logit(mu)) *
    stats::rbinom(length(mu), trials, boot::inv.logit(occ))
}

#' A posterior epred function for the rep-constant occupancy model, sufficient 
#' for \code{brms::posterior_epred()} to work. 
#' @param prep Output of \code{brms::prepare_predictions}. See brms custom 
#' families vignette at 
#' https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html
#' @return Posterior epreds
#' @noRd
posterior_epred_occupancy_single_C <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  occ <- brms::get_dpar(prep, "occ")
  trials <- prep$data$vint1
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * occ * trials
}