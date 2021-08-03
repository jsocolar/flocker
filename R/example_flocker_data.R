#' Create example data for use with \code{make_flocker_data()} and downstream functions.
#' @param visit_constant logical: return data with constant detection covariates
#' within sites (TRUE) or variable detection covariates (FALSE)
#' @return A three element named list with the observation matrix ($obs), the
#' site covariate dataframe ($site_covs), and the visit covariate list
#' ($visit_covs). If visit_constant is TRUE, then $visit_covs will be NULL.
#' @export

example_flocker_data <- function(visit_constant = FALSE, seed = 123) {
  if (!is.null(seed)) {set.seed(seed)}
  npt <- 100
  nsp <- 30
  nsite <- npt*nsp
  nvisit <- 4
  site_covs <- data.frame(sc1 = stats::rnorm(nsite), sc2 = stats::rnorm(nsite),
                          grp = sample(c(1:20), nsite, replace = T),
                          species = rep(paste0("sp_", 1:nsp), npt))
  logit_psi <- 0 + 1*site_covs$sc1 - 1*site_covs$sc2 + rnorm(20)[site_covs$grp]
  true_Z <- stats::rbinom(nsite, 1, boot::inv.logit(logit_psi))
  obs <- matrix(0, nrow=nsite, ncol=nvisit)
  if (visit_constant) {
    visit_covs = NULL
    for (i in 1:nsite) {
      if (true_Z[i] == 1) {
        obs[i, ] <- rbinom(4, 1, .5)
      }
    }
  } else {
    visit_covs <- list(vc1 = matrix(stats::rnorm(nsite*nvisit), nrow=nsite), vc2 = matrix(stats::rnorm(nsite*nvisit), nrow=nsite))
    logit_theta <- matrix(rep((-.5 + .5*site_covs$sc1 + rnorm(20)[site_covs$grp]), 4), nrow = nsite) +
      .5*visit_covs$vc1 - 1*visit_covs$vc2
    for (i in 1:nsite) {
      if (true_Z[i] == 1) {
        for (j in 1:4) {
          obs[i, j] <- stats::rbinom(1, 1, boot::inv.logit(logit_theta[i, j]))
        }
      }
    }
  }
  example_data <- list(obs = obs, site_covs = site_covs, visit_covs = visit_covs)
  example_data
}


