#' Create example data for use with \code{make_flocker_data()} and downstream functions.
#' @param rep_constant logical: create data with unit covariates only (TRUE) 
#' or data that includes event covariates (FALSE)
#' @param seed random seed. To use existing RNG state set to NULL.
#' @param n_pt number of points to simulate.
#' @param n_sp number of species to simulate.
#' @param n_rep number of replicates to simulate.
#' @return A three element named list with the observation matrix ($obs), the
#' unit covariate dataframe ($unit_covs), and the event covariate list
#' ($rep_covs). If rep_constant is TRUE, then $rep_covs will be NULL.
#' @export

example_flocker_data <- function(rep_constant = FALSE, seed = 123, 
                                 n_pt = 30, n_sp = 30, n_rep = 4) {
  if (!is.null(seed)) {set.seed(seed)}
  n_unit <- n_pt*n_sp
  
  backbone <- expand.grid(species = factor(paste0("sp_", 1:n_sp)), 
                          id_rep = 1:n_rep,
                          id_point = 1:n_pt)
  
  unit_covs <- within(unique(backbone[c("species", "id_point")]), {
    id_unit = 1:n_unit
    uc1 = stats::rnorm(n_pt)[id_point]
    uc2 = stats::rnorm(n_pt)[id_point]
    grp = sample(c(1:20), n_unit, replace = T)
    logit_psi = 0 + 1*uc1 - 1*uc2 + stats::rnorm(20)[grp]
    true_Z = stats::rbinom(n_unit, 1, boot::inv.logit(logit_psi))
  })
  
  df_full <- within(merge(backbone, unit_covs), {
    id_point_rep <- interaction(id_point, id_rep)
    ec1 = stats::rnorm(n_pt * n_rep)[id_point_rep]
    ec2 = stats::rnorm(n_pt * n_rep)[id_point_rep]
    logit_theta = -.5 + .5*uc1 + stats::rnorm(20)[grp]
    if (!rep_constant) {
      logit_theta = logit_theta + .5*ec1 - 1*ec2
    }
    obs <- stats::rbinom(n_unit * n_rep, true_Z, boot::inv.logit(logit_theta))
  })
  
  # format data for return
  out <- list(obs = t(unstack(df_full[c("obs", "id_unit")], obs ~ id_unit)), 
              unit_covs = unit_covs[c("uc1", "uc2", "grp", "species")], 
              event_covs = list(
                ec1 = t(unstack(df_full[c("ec1", "id_unit")], ec1 ~ id_unit)), 
                ec2 = t(unstack(df_full[c("ec2", "id_unit")], ec2 ~ id_unit))
              ))
  rownames(out$obs) <- NULL
  rownames(out$event_covs$ec1) <- NULL
  rownames(out$event_covs$ec2) <- NULL
  return(out)
}
