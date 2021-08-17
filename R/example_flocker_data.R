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
  
  backbone <- expand.grid(id_species = factor(paste0("sp_", 1:n_sp)), 
                          id_visit = 1:n_visit,
                          id_point = 1:n_pt)
  
  unit_covs <- within(unique(backbone[c("id_species", "id_point")]), {
    id_unit = 1:n_unit
    sc1 = stats::rnorm(n_pt)[id_point]
    sc2 = stats::rnorm(n_pt)[id_point]
    grp = sample(c(1:20), n_unit, replace = T)
    logit_psi = 0 + 1*sc1 - 1*sc2 + stats::rnorm(20)[grp]
    true_Z = stats::rbinom(n_unit, 1, boot::inv.logit(logit_psi))
  })
  
  df_full <- within(merge(backbone, unit_covs), {
    id_point_visit <- interaction(id_point, id_visit)
    vc1 = stats::rnorm(n_pt * n_visit)[id_point_visit]
    vc2 = stats::rnorm(n_pt * n_visit)[id_point_visit]
    logit_theta = -.5 + .5*sc1 + stats::rnorm(20)[grp] + .5*vc1 - 1*vc2
    obs = if(visit_constant) {
      stats::rbinom(n_unit * n_visit, true_Z, .5) 
    } else {
      stats::rbinom(n_unit * n_visit, true_Z, boot::inv.logit(logit_theta))
    }
  })
  
  # format data for return
  out <- list(obs = t(unstack(df_full[c("obs", "id_unit")], obs ~ id_unit)), 
              unit_covs = unit_covs[c("sc1", "sc2", "grp", "id_species")], 
              visit_covs = list(
                vc1 = t(unstack(df_full[c("vc1", "id_unit")], vc1 ~ id_unit)), 
                vc2 = t(unstack(df_full[c("vc2", "id_unit")], vc2 ~ id_unit))
              ))
  rownames(out$obs) <- NULL
  rownames(out$visit_covs$vc1) <- NULL
  rownames(out$visit_covs$vc2) <- NULL
  return(out)
}

