#' Create example data for use with \code{make_flocker_data()} and downstream functions.
#' @param visit_constant logical: return data with constant detection covariates
#' within units (TRUE) or variable detection covariates (FALSE)
#' @param seed random seed
#' @return A three element named list with the observation matrix ($obs), the
#' unit covariate dataframe ($unit_covs), and the visit covariate list
#' ($visit_covs). If visit_constant is TRUE, then $visit_covs will be NULL.
#' @export

example_flocker_data <- function(visit_constant = FALSE, seed = 123) {
  if (!is.null(seed)) {set.seed(seed)}
  n_pt <- 100
  n_sp <- 30
  n_unit <- n_pt*n_sp
  n_visit <- 4
  unit_covs <- data.frame(sc1 = stats::rnorm(n_unit), sc2 = stats::rnorm(n_unit),
                          grp = sample(c(1:20), n_unit, replace = T),
                          species = rep(paste0("sp_", 1:n_sp), n_pt))
  logit_psi <- 0 + 1*unit_covs$sc1 - 1*unit_covs$sc2 + stats::rnorm(20)[unit_covs$grp]
  true_Z <- stats::rbinom(n_unit, 1, boot::inv.logit(logit_psi))
  obs <- matrix(0, nrow=n_unit, ncol=n_visit)
  if (visit_constant) {
    visit_covs = NULL
    for (i in 1:n_unit) {
      if (true_Z[i] == 1) {
        obs[i, ] <- stats::rbinom(4, 1, .5)
      }
    }
  } else {
    visit_covs <- list(vc1 = matrix(stats::rnorm(n_unit*n_visit), nrow=n_unit), vc2 = matrix(stats::rnorm(n_unit*n_visit), nrow=n_unit))
    logit_theta <- matrix(rep((-.5 + .5*unit_covs$sc1 + stats::rnorm(20)[unit_covs$grp]), 4), nrow = n_unit) +
      .5*visit_covs$vc1 - 1*visit_covs$vc2
    for (i in 1:n_unit) {
      if (true_Z[i] == 1) {
        for (j in 1:4) {
          obs[i, j] <- stats::rbinom(1, 1, boot::inv.logit(logit_theta[i, j]))
        }
      }
    }
  }
  example_data <- list(obs = obs, unit_covs = unit_covs, visit_covs = visit_covs)
  example_data
}


